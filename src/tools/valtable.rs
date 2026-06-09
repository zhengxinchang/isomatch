use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::PathBuf;

use flate2::{Compression, write::GzEncoder};
use log::{info, warn};

use num_format::{Locale, ToFormattedString};
use serde::Serialize;

use crate::{
    ValTableArgs,
    tools::tools_error::ToolError,
    traits::ArgValidate,
    utils::{greetings2, open_file_bufread, print_json_block, require_file},
};

impl ArgValidate for ValTableArgs {
    fn validate(&self) {
        let mut error_msg = String::new();
        let mut has_error = false;

        require_file(
            "Merged GTF",
            &self.query_gtf,
            &mut error_msg,
            &mut has_error,
        );
        for input in &self.inputs {
            require_file("Source GTF", input, &mut error_msg, &mut has_error);
        }

        if has_error {
            panic!("{}", error_msg);
        }
    }
}

#[derive(Serialize)]
struct SampleRecord {
    file_id: u32,
    file_name: String,
}

#[derive(Serialize)]
struct FileStats {
    file_name: String,
    extracted: u32,
    attr_missing: u32,
    src_tx_not_in_merged: u32,
}

#[derive(Serialize)]
struct ValtableStats {
    merged_tx_count: usize,
    merged_gtf_sample_count: usize,
    merged_gtf_samples: Vec<SampleRecord>,
    per_file: Vec<FileStats>,
}

pub fn run_valtable(args: &ValTableArgs) -> Result<(), ToolError> {
    greetings2(args);
    args.validate();

    let default_val = args.default_val.as_deref().unwrap_or("0.0");

    let mut sample_map: HashMap<u32, String> = HashMap::new();
    let mut projection: HashMap<(String, u32), usize> = HashMap::new();
    let mut isomatch_txids: Vec<String> = Vec::new();
    let mut is_isom_gtf = false;

    info!("Loading merged GTF: {}", args.query_gtf.display());
    let mut reader = open_file_bufread(&args.query_gtf)?;
    let mut line = String::new();
    let mut line_no = 0;

    while reader.read_line(&mut line)? != 0 {
        line_no += 1;
        if line.starts_with("##ISOM") {
            is_isom_gtf = true;
            if line.starts_with("##ISOM <SAMPLE>") {
                // ##ISOM <SAMPLE> id="S1"; input="path/to/file.gtf.gz";
                let id_str = parse_kv_quoted(&line, "id");
                let input_path = parse_kv_quoted(&line, "input");
                match (id_str, input_path) {
                    (None, _) => warn!(
                        "ISOM <SAMPLE> line missing 'id' attribute, skipping: {}",
                        line.trim_end()
                    ),
                    (_, None) => warn!(
                        "ISOM <SAMPLE> line missing 'input' attribute, skipping: {}",
                        line.trim_end()
                    ),
                    (Some(id_str), Some(input_path)) => {
                        match id_str.trim_start_matches('S').parse::<u32>() {
                            Err(_) => warn!(
                                "ISOM <SAMPLE> line has unparseable file_id '{}', skipping: {}",
                                id_str,
                                line.trim_end()
                            ),
                            Ok(0) => warn!(
                                "ISOM <SAMPLE> line has invalid file_id 0, skipping: {}",
                                line.trim_end()
                            ),
                            Ok(file_id) => {
                                let base_name = PathBuf::from(&input_path)
                                    .file_name()
                                    .map(|f| f.to_string_lossy().to_string())
                                    .unwrap_or(input_path);
                                sample_map.insert(file_id, base_name);
                            }
                        }
                    }
                }
            }
            line.clear();
            continue;
        }

        if line.starts_with('#') {
            line.clear();
            continue;
        }

        let cols: Vec<&str> = line.splitn(9, '\t').collect();
        if cols.len() < 9 {
            warn!("Bad record line at {}", line_no);
            line.clear();
            continue;
        }

        if cols[2] != "transcript" {
            line.clear();
            continue;
        }

        let attrs = cols[8];
        let tx_id = extract_gtf_attr(attrs, "transcript_id");
        if tx_id.is_empty() {
            return Err(ToolError::ReadMergedGTFFailed {
                reason: format!("Can not find transcript id in line {}", line_no),
            });
        }

        let tx_index = isomatch_txids.len();
        isomatch_txids.push(tx_id);

        // ISOM_SRC: "S1:src_tx_id:start:end:...|S2:src_tx_id:..."
        let isom_src = extract_gtf_attr(attrs, "ISOM_SRC");
        for src_record in isom_src.split('|') {
            if src_record.is_empty() {
                return Err(ToolError::ReadMergedGTFFailed {
                    reason: format!("Can not read ISOM_SRC from line: {}", line_no),
                });
            }
            let mut parts = src_record.splitn(3, ':');
            let file_id_str = parts.next().unwrap_or("");
            let src_tx_id = parts.next().unwrap_or("").to_string();
            if file_id_str.is_empty() || src_tx_id.is_empty() {
                return Err(ToolError::ReadMergedGTFFailed {
                    reason: format!(
                        "Can not read file id string from ISOM_SRC in line: {}",
                        line_no
                    ),
                });
            }
            let file_id: u32 = match file_id_str.trim_start_matches('S').parse::<u32>() {
                Ok(0) => {
                    return Err(ToolError::ReadMergedGTFFailed {
                        reason: format!("ISOM_SRC has invalid file_id 0 in line {}", line_no),
                    });
                }
                Ok(id) => id,
                Err(_) => {
                    return Err(ToolError::ReadMergedGTFFailed {
                        reason: format!(
                            "ISOM_SRC has unparseable file_id '{}' in line {}",
                            file_id_str, line_no
                        ),
                    });
                }
            };
            projection.insert((src_tx_id, file_id), tx_index);
        }

        line.clear();
    }

    if !is_isom_gtf {
        return Err(ToolError::ReadMergedGTFFailed {
            reason: format!(
                "{} has no ISOM headers; is this an isomatch merged GTF?",
                args.query_gtf.display()
            ),
        });
    }

    let n_txs = isomatch_txids.len();
    info!(
        "Loaded {} merged transcripts, {} sample(s) from ISOM headers",
        n_txs.to_formatted_string(&Locale::en),
        sample_map.len().to_formatted_string(&Locale::en)
    );

    // Process each source GTF: extract attr_val, write one column to a temp file.
    let tmp_dir = std::env::temp_dir();
    let out_stem = args
        .out
        .file_name()
        .map(|f| f.to_string_lossy().to_string())
        .unwrap_or_else(|| "valtable".to_string());

    let mut temp_file_paths: Vec<PathBuf> = Vec::new();
    let mut column_headers: Vec<String> = Vec::new();
    let mut all_stats: Vec<FileStats> = Vec::new();

    for src_path in &args.inputs {
        let base_name = src_path
            .file_name()
            .map(|f| f.to_string_lossy().to_string())
            .ok_or_else(|| ToolError::InvalidPath {
                path: src_path.to_string_lossy().to_string(),
            })?;

        // Match the source file's base name against ISOM header entries to get file_id.
        let file_id = sample_map
            .iter()
            .find(|(_, v)| v.as_str() == base_name)
            .map(|(k, _)| *k);

        let Some(file_id) = file_id else {
            warn!(
                "Source GTF '{}' not found in merged GTF ISOM headers, skipping",
                base_name
            );
            continue;
        };

        info!("Processing S{}: {}", file_id, base_name);

        let mut col: Vec<String> = vec![default_val.to_string(); n_txs];
        let mut stats = FileStats {
            file_name: base_name.clone(),
            extracted: 0,
            attr_missing: 0,
            src_tx_not_in_merged: 0,
        };

        let mut src_reader = open_file_bufread(src_path)?;
        let mut src_line = String::new();
        let mut src_line_no = 0u32;
        while src_reader.read_line(&mut src_line)? != 0 {
            src_line_no += 1;
            if src_line.starts_with('#') {
                src_line.clear();
                continue;
            }
            let cols: Vec<&str> = src_line.splitn(9, '\t').collect();
            if cols.len() < 9 {
                warn!(
                    "{} line {}: fewer than 9 columns, skipping",
                    base_name, src_line_no
                );
                src_line.clear();
                continue;
            }
            if cols[2] != "transcript" {
                src_line.clear();
                continue;
            }

            let attrs = cols[8];
            let src_tx_id = extract_gtf_attr(attrs, "transcript_id");
            if src_tx_id.is_empty() {
                warn!(
                    "{} line {}: transcript line missing transcript_id, skipping",
                    base_name, src_line_no
                );
                src_line.clear();
                continue;
            }

            if let Some(&tx_index) = projection.get(&(src_tx_id, file_id)) {
                let attr_val = extract_gtf_attr(attrs, &args.attr_val);
                if attr_val.is_empty() {
                    stats.attr_missing += 1;
                } else {
                    col[tx_index] = attr_val;
                    stats.extracted += 1;
                }
            } else {
                stats.src_tx_not_in_merged += 1;
            }

            src_line.clear();
        }

        let tmp_path = tmp_dir.join(format!("{}_{}.tmp", out_stem, file_id));
        let mut tmp_writer = BufWriter::new(File::create(&tmp_path)?);
        for val in &col {
            tmp_writer.write_all(val.as_bytes())?;
            tmp_writer.write_all(b"\n")?;
        }
        tmp_writer.flush()?;

        temp_file_paths.push(tmp_path);
        column_headers.push(base_name);
        all_stats.push(stats);
    }

    // Combine temp files into final matrix (transcript_id + one column per source file).
    let mut out_path = args.out.clone();
    out_path.add_extension("valtable.tsv.gz");
    let mut out_writer = BufWriter::new(GzEncoder::new(
        File::create(&out_path)?,
        Compression::default(),
    ));

    write!(out_writer, "transcript_id")?;
    for header in &column_headers {
        write!(out_writer, "\t{}", header)?;
    }
    writeln!(out_writer)?;

    let columns: Vec<Vec<String>> = temp_file_paths
        .iter()
        .map(|p| {
            let rdr = open_file_bufread(p)?;
            rdr.lines().collect::<std::io::Result<_>>()
        })
        .collect::<std::io::Result<_>>()?;

    for (i, tx_id) in isomatch_txids.iter().enumerate() {
        write!(out_writer, "{}", tx_id)?;
        for col in &columns {
            let val = col.get(i).map(|s| s.as_str()).unwrap_or(default_val);
            write!(out_writer, "\t{}", val)?;
        }
        writeln!(out_writer)?;
    }
    out_writer.flush()?;

    for p in &temp_file_paths {
        let _ = std::fs::remove_file(p);
    }

    let mut merged_gtf_samples: Vec<SampleRecord> = sample_map
        .iter()
        .map(|(&file_id, file_name)| SampleRecord {
            file_id,
            file_name: file_name.clone(),
        })
        .collect();
    merged_gtf_samples.sort_by_key(|s| s.file_id);

    let stats = ValtableStats {
        merged_tx_count: n_txs,
        merged_gtf_sample_count: merged_gtf_samples.len(),
        merged_gtf_samples,
        per_file: all_stats,
    };

    print_json_block("Valtable stats", &stats);

    let mut stats_path = args.out.clone();
    stats_path.add_extension("valtable_stats.json");
    let stats_json = serde_json::to_string_pretty(&stats)?;
    std::fs::write(&stats_path, stats_json)?;

    info!("Output saved to: {}", out_path.display());
    info!("Stats saved to: {}", stats_path.display());
    info!("Finished!");
    Ok(())
}

// Parse a key="value" pair from an ISOM header line.
fn parse_kv_quoted(line: &str, key: &str) -> Option<String> {
    let needle = format!("{}=\"", key);
    let start = line.find(&needle)? + needle.len();
    let end = line[start..].find('"')? + start;
    Some(line[start..end].to_string())
}

// Extract the value of a named attribute from a GTF attribute column string.
fn extract_gtf_attr(attrs: &str, key: &str) -> String {
    for attr in attrs.split(';') {
        let attr = attr.trim();
        if attr.is_empty() {
            continue;
        }
        let key_part = attr.split_ascii_whitespace().next().unwrap_or("");
        if key_part != key {
            continue;
        }
        if let Some(q_start) = attr.find('"') {
            if let Some(q_len) = attr[q_start + 1..].find('"') {
                return attr[q_start + 1..q_start + 1 + q_len].to_string();
            }
        }
        if let Some(val) = attr.split_ascii_whitespace().nth(1) {
            return val.to_string();
        }
    }
    String::new()
}
