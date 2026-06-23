//! revert isomatch merged sample into multiple GTFs

use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufWriter, Write},
    path::{Path, PathBuf},
};

use flate2::{Compression, write::GzEncoder};
use log::{info, warn};
use serde::Serialize;

use crate::{
    RevertArgs,
    index::gtf::parse_gtf_attr_value,
    tools::tools_error::ToolError,
    traits::ArgValidate,
    utils::{
        IsomSample, greetings2, is_gzipped, is_isomatch_merged_gtf, open_file_bufread,
        print_json_block, read_samples, require_file,
    },
};

const MAX_OPEN_FILES: usize = 512;

impl ArgValidate for RevertArgs {
    fn validate(&self) {
        let mut error_msg = String::new();
        let mut has_error = false;

        require_file("Input GTF", &self.input, &mut error_msg, &mut has_error);

        if let Some(track_file) = &self.track_file {
            require_file("Track TSV", track_file, &mut error_msg, &mut has_error);
        }

        let out_parent = self.out.parent().unwrap_or_else(|| Path::new("."));
        if !out_parent.exists() {
            error_msg.push_str(&format!(
                "\nOutput directory parent does not exist: {:?}",
                out_parent
            ));
            has_error = true;
        } else if !out_parent.is_dir() {
            error_msg.push_str(&format!(
                "\nOutput directory parent is not a directory: {:?}",
                out_parent
            ));
            has_error = true;
        }

        match is_isomatch_merged_gtf(&self.input) {
            Ok(true) => {}
            Ok(false) => {
                error_msg.push_str(&format!(
                    "\nInput has no matching ISOM schema header: {:?}",
                    self.input
                ));
                has_error = true;
            }
            Err(err) => {
                error_msg.push_str(&format!(
                    "\nCan not inspect input GTF header {:?}: {}",
                    self.input, err
                ));
                has_error = true;
            }
        }

        if has_error {
            panic!("{}", error_msg);
        }
    }
}

pub fn run_revert(args: &RevertArgs) -> Result<(), ToolError> {
    greetings2(&args);
    args.validate();

    let samples = read_samples(&args.input)?;
    if samples.is_empty() {
        return Err(ToolError::ReadMergedGTFFailed {
            reason: "No ##ISOM <SAMPLE> header found".to_string(),
        });
    }

    std::fs::create_dir_all(&args.out)?;

    let track_gene_ids = match &args.track_file {
        Some(track_file) => Some(read_track_gene_ids(track_file, &samples)?),
        None => None,
    };

    let mut file_ids: Vec<u32> = samples.keys().copied().collect();
    file_ids.sort_unstable();

    let mut stats = RevertStats::new(args, &samples);

    for (batch_idx, chunk) in file_ids.chunks(MAX_OPEN_FILES).enumerate() {
        info!("Reverting sample batch with {} output file(s)", chunk.len());
        let mut writers = open_writers(&samples, chunk, args)?;
        write_reverted_batch(
            &args.input,
            &mut writers,
            track_gene_ids.as_ref(),
            batch_idx == 0,
            &mut stats,
        )?;
        for writer in writers.values_mut() {
            writer.flush()?;
        }
    }

    print_json_block("Revert stats", &stats);
    let stats_path = args.out.join("revert_stats.json");
    let stats_json = serde_json::to_string_pretty(&stats)?;
    std::fs::write(&stats_path, stats_json)?;

    info!("Reverted GTFs saved to: {}", args.out.display());
    info!("Stats saved to: {}", stats_path.display());
    info!("Finished!");
    Ok(())
}

#[derive(Serialize)]
struct RevertFileStats {
    file_id: u32,
    file_name: String,
    output_path: String,
    transcripts: u64,
    exons: u64,
    missing_track_gene_id: u64,
}

#[derive(Serialize)]
struct RevertStats {
    merged_gtf: String,
    output_dir: String,
    track_file: Option<String>,
    gzipped_output: bool,
    merged_tx_count: u64,
    source_tx_count: u64,
    source_exon_count: u64,
    output_gtf_count: usize,
    missing_track_gene_id: u64,
    per_file: Vec<RevertFileStats>,

    #[serde(skip)]
    per_file_idx: HashMap<u32, usize>,
}

impl RevertStats {
    fn new(args: &RevertArgs, samples: &HashMap<u32, IsomSample>) -> Self {
        let mut file_ids: Vec<u32> = samples.keys().copied().collect();
        file_ids.sort_unstable();

        let mut per_file_idx = HashMap::new();
        let per_file = file_ids
            .iter()
            .enumerate()
            .map(|(idx, &file_id)| {
                per_file_idx.insert(file_id, idx);
                let sample = &samples[&file_id];
                RevertFileStats {
                    file_id,
                    file_name: sample.file_name.clone(),
                    output_path: output_path_for_sample(&args.out, sample, args)
                        .display()
                        .to_string(),
                    transcripts: 0,
                    exons: 0,
                    missing_track_gene_id: 0,
                }
            })
            .collect();

        Self {
            merged_gtf: args.input.display().to_string(),
            output_dir: args.out.display().to_string(),
            track_file: args
                .track_file
                .as_ref()
                .map(|track_file| track_file.display().to_string()),
            gzipped_output: args.gzipped,
            merged_tx_count: 0,
            source_tx_count: 0,
            source_exon_count: 0,
            output_gtf_count: samples.len(),
            missing_track_gene_id: 0,
            per_file,
            per_file_idx,
        }
    }

    fn observe_merged_tx(&mut self) {
        self.merged_tx_count += 1;
    }

    fn observe_source_tx(&mut self, file_id: u32, exon_count: usize) {
        self.source_tx_count += 1;
        self.source_exon_count += exon_count as u64;
        if let Some(&idx) = self.per_file_idx.get(&file_id) {
            self.per_file[idx].transcripts += 1;
            self.per_file[idx].exons += exon_count as u64;
        }
    }

    fn observe_missing_track_gene_id(&mut self, file_id: u32) {
        self.missing_track_gene_id += 1;
        if let Some(&idx) = self.per_file_idx.get(&file_id) {
            self.per_file[idx].missing_track_gene_id += 1;
        }
    }
}

#[derive(Debug)]
struct MergedBlock {
    chrom: String,
    strand: String,
    merged_tx_id: String,
    merged_gene_id: String,
    sources: Vec<SourceRecord>,
    exons: Vec<(u32, u32)>,
}

impl MergedBlock {
    fn from_gtf_cols(cols: &[&str], line_no: usize) -> Result<Self, ToolError> {
        let attrs = cols[8];
        let merged_tx_id = require_attr(attrs, "transcript_id", line_no)?;
        let merged_gene_id = require_attr(attrs, "gene_id", line_no)?;
        let isom_src = require_attr(attrs, "ISOM_SRC", line_no)?;

        Ok(Self {
            chrom: cols[0].to_string(),
            strand: cols[6].to_string(),
            merged_tx_id,
            merged_gene_id,
            sources: parse_sources(&isom_src, line_no)?,
            exons: Vec::new(),
        })
    }

    fn write_sources(
        &self,
        writers: &mut HashMap<u32, Box<dyn Write>>,
        track_gene_ids: Option<&TrackGeneIds>,
        stats: &mut RevertStats,
    ) -> Result<(), ToolError> {
        if self.exons.is_empty() {
            return Err(ToolError::ReadMergedGTFFailed {
                reason: format!("Merged transcript {} has no exon", self.merged_tx_id),
            });
        }

        for source in &self.sources {
            let Some(writer) = writers.get_mut(&source.file_id) else {
                continue;
            };

            let gene_id = match track_gene_ids {
                Some(gene_ids) => gene_ids
                    .get(&(
                        source.file_id,
                        self.merged_tx_id.clone(),
                        source.tx_id.clone(),
                    ))
                    .map(String::as_str)
                    .unwrap_or_else(|| {
                        warn!(
                            "No track gene_id for S{}:{} in {}, using merged gene_id",
                            source.file_id, source.tx_id, self.merged_tx_id
                        );
                        stats.observe_missing_track_gene_id(source.file_id);
                        self.merged_gene_id.as_str()
                    }),
                None => self.merged_gene_id.as_str(),
            };

            let source_exons = source.exons_from_repr(&self.exons)?;
            stats.observe_source_tx(source.file_id, source_exons.len());

            write_gtf_record(
                writer.as_mut(),
                &self.chrom,
                "transcript",
                source.start,
                source.end,
                &self.strand,
                gene_id,
                &source.tx_id,
                &self.merged_tx_id,
                None,
            )?;

            for (idx, (start, end)) in source_exons.iter().enumerate() {
                write_gtf_record(
                    writer.as_mut(),
                    &self.chrom,
                    "exon",
                    *start,
                    *end,
                    &self.strand,
                    gene_id,
                    &source.tx_id,
                    &self.merged_tx_id,
                    Some(idx + 1),
                )?;
            }
        }

        Ok(())
    }
}

#[derive(Debug)]
struct SourceRecord {
    file_id: u32,
    tx_id: String,
    start: u32,
    end: u32,
    junction_diffs: String,
}

impl SourceRecord {
    fn exons_from_repr(&self, repr_exons: &[(u32, u32)]) -> Result<Vec<(u32, u32)>, ToolError> {
        let mut exons = repr_exons.to_vec();
        exons[0].0 = self.start;
        let last_idx = exons.len() - 1;
        exons[last_idx].1 = self.end;

        for (junction_no, left_diff, right_diff) in parse_junction_diffs(&self.junction_diffs)? {
            if junction_no == 0 || junction_no >= exons.len() {
                return Err(ToolError::ReadMergedGTFFailed {
                    reason: format!(
                        "Invalid junction diff index {junction_no} for {}",
                        self.tx_id
                    ),
                });
            }
            let junction_idx = junction_no - 1;
            exons[junction_idx].1 = shift_coordinate(repr_exons[junction_idx].1, left_diff)?;
            exons[junction_idx + 1].0 =
                shift_coordinate(repr_exons[junction_idx + 1].0, right_diff)?;
        }

        Ok(exons)
    }
}

type TrackGeneIds = HashMap<(u32, String, String), String>;

fn read_track_gene_ids(
    track_file: &Path,
    samples: &HashMap<u32, IsomSample>,
) -> Result<TrackGeneIds, ToolError> {
    let sample_file_ids: HashMap<&str, u32> = samples
        .iter()
        .map(|(&file_id, sample)| (sample.file_name.as_str(), file_id))
        .collect();

    let mut reader = open_file_bufread(track_file)?;
    let mut line = String::new();
    if reader.read_line(&mut line)? == 0 {
        return Err(ToolError::ReadMergedGTFFailed {
            reason: format!("Track file is empty: {}", track_file.display()),
        });
    }

    let header: Vec<&str> = line.trim_end_matches(['\r', '\n']).split('\t').collect();
    let merged_tx_idx = track_col(&header, "merged_tx_id")?;
    let src_tx_idx = track_col(&header, "src_tx_id")?;
    let src_gene_idx = track_col(&header, "src_gene_id")?;
    let src_file_idx = track_col(&header, "src_file_name")?;
    let max_idx = merged_tx_idx
        .max(src_tx_idx)
        .max(src_gene_idx)
        .max(src_file_idx);

    let mut out = HashMap::new();
    line.clear();
    while reader.read_line(&mut line)? != 0 {
        let cols: Vec<&str> = line.trim_end_matches(['\r', '\n']).split('\t').collect();
        if cols.len() <= max_idx {
            line.clear();
            continue;
        }

        if let Some(&file_id) = sample_file_ids.get(cols[src_file_idx]) {
            out.insert(
                (
                    file_id,
                    cols[merged_tx_idx].to_string(),
                    cols[src_tx_idx].to_string(),
                ),
                cols[src_gene_idx].to_string(),
            );
        }
        line.clear();
    }

    Ok(out)
}

fn track_col(header: &[&str], name: &str) -> Result<usize, ToolError> {
    header
        .iter()
        .position(|col| *col == name)
        .ok_or_else(|| ToolError::ReadMergedGTFFailed {
            reason: format!("Track file missing column: {name}"),
        })
}

fn open_writers(
    samples: &HashMap<u32, IsomSample>,
    file_ids: &[u32],
    args: &RevertArgs,
) -> Result<HashMap<u32, Box<dyn Write>>, ToolError> {
    let mut writers = HashMap::new();
    for file_id in file_ids {
        let sample = samples
            .get(file_id)
            .ok_or_else(|| ToolError::ReadMergedGTFFailed {
                reason: format!("No sample header found for S{file_id}"),
            })?;

        let out_path = output_path_for_sample(&args.out, sample, args);
        let file = File::create(&out_path)?;
        let writer: Box<dyn Write> = if is_gzipped(&out_path) {
            Box::new(BufWriter::new(GzEncoder::new(file, Compression::default())))
        } else {
            Box::new(BufWriter::new(file))
        };
        writers.insert(*file_id, writer);
    }
    Ok(writers)
}

fn output_path_for_sample(out_dir: &Path, sample: &IsomSample, args: &RevertArgs) -> PathBuf {
    let mut out_path = out_dir.join(&sample.file_name);
    if args.gzipped && !is_gzipped(&out_path) {
        out_path.add_extension("gz");
    }
    out_path
}

/// take input merged gtf and a set of activated file id to revert them in to single files.
fn write_reverted_batch(
    merged_gtf: &Path,
    writers: &mut HashMap<u32, Box<dyn Write>>,
    track_gene_ids: Option<&TrackGeneIds>,
    observe_merged_txs: bool,
    stats: &mut RevertStats,
) -> Result<(), ToolError> {
    let mut reader = open_file_bufread(merged_gtf)?;
    let mut line = String::new();
    let mut block: Option<MergedBlock> = None;
    let mut line_no = 0usize;

    while reader.read_line(&mut line)? != 0 {
        line_no += 1;
        if line.starts_with('#') {
            line.clear();
            continue;
        }

        let line_trimmed = line.trim_end_matches(['\r', '\n']);
        let cols: Vec<&str> = line_trimmed.splitn(9, '\t').collect();
        if cols.len() < 9 {
            return Err(ToolError::ReadMergedGTFFailed {
                reason: format!("Line {line_no} has fewer than 9 columns"),
            });
        }

        match cols[2] {
            "transcript" => {
                if let Some(prev_block) = block.take() {
                    if observe_merged_txs {
                        stats.observe_merged_tx();
                    }
                    prev_block.write_sources(writers, track_gene_ids, stats)?;
                }
                block = Some(MergedBlock::from_gtf_cols(&cols, line_no)?);
            }
            "exon" => {
                let Some(curr_block) = block.as_mut() else {
                    // because isomatch merged GTF must have a sorted transcript
                    // block. if not, then must have error here.
                    return Err(ToolError::ReadMergedGTFFailed {
                        reason: format!("Exon before transcript at line {line_no}"),
                    });
                };
                curr_block
                    .exons
                    .push((cols[3].parse::<u32>()?, cols[4].parse::<u32>()?));
            }
            _ => {}
        }

        line.clear();
    }

    // process the last available block
    if let Some(prev_block) = block.take() {
        if observe_merged_txs {
            stats.observe_merged_tx();
        }
        prev_block.write_sources(writers, track_gene_ids, stats)?;
    }

    Ok(())
}

fn parse_sources(isom_src: &str, line_no: usize) -> Result<Vec<SourceRecord>, ToolError> {
    let mut out = Vec::new();
    for record in isom_src.split('|') {
        let parts: Vec<&str> = record.splitn(8, ':').collect();
        if parts.len() != 8 {
            return Err(ToolError::ReadMergedGTFFailed {
                reason: format!("Malformed ISOM_SRC record at line {line_no}: {record}"),
            });
        }
        out.push(SourceRecord {
            file_id: parse_source_file_id(parts[0])?,
            tx_id: parts[1].to_string(),
            start: parts[2].parse::<u32>()?,
            end: parts[3].parse::<u32>()?,
            junction_diffs: parts[7].to_string(),
        });
    }
    Ok(out)
}

fn parse_junction_diffs(junction_diffs: &str) -> Result<Vec<(usize, i32, i32)>, ToolError> {
    if junction_diffs == "no_diff" {
        return Ok(Vec::new());
    }

    let mut out = Vec::new();
    let mut rest = junction_diffs.trim();
    while !rest.is_empty() {
        let start = rest
            .find('(')
            .ok_or_else(|| ToolError::ReadMergedGTFFailed {
                reason: format!("Malformed junction diff list: {junction_diffs}"),
            })?;
        let end = rest[start..]
            .find(')')
            .map(|idx| start + idx)
            .ok_or_else(|| ToolError::ReadMergedGTFFailed {
                reason: format!("Malformed junction diff list: {junction_diffs}"),
            })?;
        out.push(parse_junction_diff(&rest[start..=end])?);
        rest = rest[end + 1..].trim_start_matches([',', ' ', '\t']);
    }

    Ok(out)
}

fn parse_junction_diff(diff: &str) -> Result<(usize, i32, i32), ToolError> {
    let diff = diff.trim();
    let diff = diff
        .strip_prefix('(')
        .and_then(|value| value.strip_suffix(')'))
        .ok_or_else(|| ToolError::ReadMergedGTFFailed {
            reason: format!("Malformed junction diff: {diff}"),
        })?;
    let parts: Vec<&str> = diff.split(',').collect();
    if parts.len() != 3 {
        return Err(ToolError::ReadMergedGTFFailed {
            reason: format!("Malformed junction diff: {diff}"),
        });
    }
    Ok((
        parts[0]
            .parse::<usize>()
            .map_err(|err| ToolError::ReadMergedGTFFailed {
                reason: format!("Can not parse junction number in diff {diff}: {err}"),
            })?,
        parts[1]
            .parse::<i32>()
            .map_err(|err| ToolError::ReadMergedGTFFailed {
                reason: format!("Can not parse left diff in {diff}: {err}"),
            })?,
        parts[2]
            .parse::<i32>()
            .map_err(|err| ToolError::ReadMergedGTFFailed {
                reason: format!("Can not parse right diff in {diff}: {err}"),
            })?,
    ))
}

fn shift_coordinate(repr_coord: u32, diff: i32) -> Result<u32, ToolError> {
    let shifted = i64::from(repr_coord) - i64::from(diff);
    u32::try_from(shifted).map_err(|_| ToolError::ReadMergedGTFFailed {
        reason: format!("Invalid reconstructed coordinate: {repr_coord} - ({diff})"),
    })
}

fn write_gtf_record(
    writer: &mut dyn Write,
    chrom: &str,
    feature: &str,
    start: u32,
    end: u32,
    strand: &str,
    gene_id: &str,
    tx_id: &str,
    isom_tx_id: &str,
    exon_number: Option<usize>,
) -> Result<(), ToolError> {
    write!(
        writer,
        "{chrom}\tisomatch\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t"
    )?;
    write_attr(writer, "gene_id", gene_id)?;
    writer.write_all(b"; ")?;
    write_attr(writer, "transcript_id", tx_id)?;
    writer.write_all(b"; ")?;
    write_attr(writer, "ISOM_TX_ID", isom_tx_id)?;
    if let Some(exon_number) = exon_number {
        writer.write_all(b"; ")?;
        write_attr(writer, "exon_number", &exon_number.to_string())?;
    }
    writer.write_all(b";\n")?;
    Ok(())
}

fn write_attr(writer: &mut dyn Write, key: &str, value: &str) -> Result<(), ToolError> {
    write!(writer, "{key} \"")?;
    for byte in value.bytes() {
        match byte {
            b'\\' => writer.write_all(b"\\\\")?,
            b'"' => writer.write_all(b"\\\"")?,
            _ => writer.write_all(&[byte])?,
        }
    }
    writer.write_all(b"\"")?;
    Ok(())
}

fn require_attr(attrs: &str, key: &str, line_no: usize) -> Result<String, ToolError> {
    parse_gtf_attr_value(attrs, key).ok_or_else(|| ToolError::ReadMergedGTFFailed {
        reason: format!("Can not find {key} in line {line_no}"),
    })
}

fn parse_source_file_id(id: &str) -> Result<u32, ToolError> {
    match id.trim_start_matches('S').parse::<u32>() {
        Ok(0) => Err(ToolError::ReadMergedGTFFailed {
            reason: "Source file id 0 is invalid".to_string(),
        }),
        Ok(file_id) => Ok(file_id),
        Err(err) => Err(ToolError::ReadMergedGTFFailed {
            reason: format!("Can not parse source file id {id}: {err}"),
        }),
    }
}
