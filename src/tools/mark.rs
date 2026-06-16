//! mark transcirpts based on critera

use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{BufRead, BufWriter, Write},
};

use flate2::{Compression, write::GzEncoder};
use log::info;
use serde::Serialize;

use crate::{
    MarkArgs,
    core::tx_strand::ISOMSTRAND,
    index::gtf::parse_gtf_attr_value,
    region::{RegionDb, RegionType},
    tools::tools_error::ToolError,
    traits::ArgValidate,
    utils::{greetings2, open_file_bufread, print_json_block, require_file},
};

impl ArgValidate for MarkArgs {
    fn validate(&self) {
        let mut error_msg = String::new();
        let mut has_error = false;

        require_file("Input GTF", &self.input, &mut error_msg, &mut has_error);
        require_file(
            "Reference GTF",
            &self.track_file,
            &mut error_msg,
            &mut has_error,
        );

        if has_error {
            log::error!("Error validating arguments: {}", error_msg);
            std::process::exit(1);
        }
    }
}

pub fn run_mark(args: &MarkArgs) -> Result<(), ToolError> {
    greetings2(&args);
    args.validate();

    let chrmap = None;
    let gene_db = RegionDb::from_gtf_gene(&args.track_file, RegionType::Gene, &chrmap)?;
    let transcripts = read_transcript_regions(args)?;

    let mut out_path = args.out.clone();
    out_path.add_extension("mark_dup_gene.tsv.gz");
    let out_file = File::create(&out_path)?;
    let mut writer = BufWriter::new(GzEncoder::new(out_file, Compression::default()));
    writeln!(
        writer,
        "transcript_id\tgene_id\tn_exons\toverlapped_gene_count\toverlapped_gene_strand_type_count\toverlapped_gene_strands\toverlapped_genes"
    )?;

    let mut overlapped_tx_count = 0usize;

    for tx in transcripts.values() {
        let hits = gene_db.query_overlaps_range_all_strands(&tx.chrom, tx.start, tx.end);
        let mut gene_hits = BTreeSet::new();
        let mut strands = BTreeSet::new();
        for hit in hits {
            strands.insert(strand_label(hit.strand));
            gene_hits.insert(format!(
                "{}:{}:{}",
                hit.id,
                hit.name,
                strand_label(hit.strand)
            ));
        }
        if !gene_hits.is_empty() {
            overlapped_tx_count += 1;
        }
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            tx.tx_id,
            tx.gene_id,
            tx.n_exons,
            gene_hits.len(),
            strands.len(),
            join_strs_or_na(strands),
            join_strings_or_na(gene_hits),
        )?;
    }
    writer.flush()?;

    let stats = MarkStats {
        transcript_count: transcripts.len(),
        overlapped_transcript_count: overlapped_tx_count,
    };
    print_json_block("Mark stats", &stats);

    let mut stats_path = args.out.clone();
    stats_path.add_extension("mark_stats.json");
    std::fs::write(&stats_path, serde_json::to_string_pretty(&stats)?)?;

    info!("Output saved to: {}", out_path.display());
    info!("Stats saved to: {}", stats_path.display());
    info!("Finished!");

    Ok(())
}

#[derive(Debug)]
struct TranscriptRegion {
    tx_id: String,
    gene_id: String,
    chrom: String,
    start: u32,
    end: u32,
    n_exons: u32,
}

#[derive(Debug, Serialize)]
struct MarkStats {
    transcript_count: usize,
    overlapped_transcript_count: usize,
}

fn read_transcript_regions(
    args: &MarkArgs,
) -> Result<BTreeMap<String, TranscriptRegion>, ToolError> {
    let mut reader = open_file_bufread(&args.input)?;
    let mut transcripts = BTreeMap::new();
    let mut line = String::new();
    let mut line_no = 0usize;

    while reader.read_line(&mut line)? != 0 {
        line_no += 1;
        let raw = line.trim_end();
        if raw.is_empty() || raw.starts_with('#') {
            line.clear();
            continue;
        }

        let fields: Vec<&str> = raw.split('\t').collect();
        if fields.len() < 9 {
            return Err(ToolError::FailedParseGTF {
                reason: format!("line {line_no}: expected 9 columns, got {}", fields.len()),
            });
        }
        if fields[2] != "exon" {
            line.clear();
            continue;
        }

        let start = fields[3].parse::<u32>()?;
        let end = fields[4].parse::<u32>()?;
        let tx_id = parse_gtf_attr_value(fields[8], "transcript_id").ok_or_else(|| {
            ToolError::FailedParseGTF {
                reason: format!("line {line_no}: missing transcript_id"),
            }
        })?;
        let gene_id = parse_gtf_attr_value(fields[8], "gene_id").unwrap_or_default();

        transcripts
            .entry(tx_id.clone())
            .and_modify(|tx: &mut TranscriptRegion| {
                tx.start = tx.start.min(start);
                tx.end = tx.end.max(end);
                tx.n_exons += 1;
                if tx.gene_id.is_empty() {
                    tx.gene_id = gene_id.clone();
                }
            })
            .or_insert_with(|| TranscriptRegion {
                tx_id,
                gene_id,
                chrom: fields[0].to_string(),
                start,
                end,
                n_exons: 1,
            });
        line.clear();
    }

    Ok(transcripts)
}

fn strand_label(strand: ISOMSTRAND) -> &'static str {
    match strand {
        ISOMSTRAND::Plus => "+",
        ISOMSTRAND::Minus => "-",
        ISOMSTRAND::Unknown => "Unk",
    }
}

fn join_strs_or_na(values: BTreeSet<&str>) -> String {
    if values.is_empty() {
        "NA".to_string()
    } else {
        values.into_iter().collect::<Vec<_>>().join(",")
    }
}

fn join_strings_or_na(values: BTreeSet<String>) -> String {
    if values.is_empty() {
        "NA".to_string()
    } else {
        values.into_iter().collect::<Vec<_>>().join(",")
    }
}
