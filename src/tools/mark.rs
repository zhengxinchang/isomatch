//! mark transcirpts based on critera

use std::{
    collections::BTreeSet,
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

    let mut tsv_path = args.out.clone();
    tsv_path.add_extension("mark_dup_gene.tsv.gz");
    let tsv_file = File::create(&tsv_path)?;
    let mut tsv_writer = BufWriter::new(GzEncoder::new(tsv_file, Compression::default()));
    writeln!(
        tsv_writer,
        "transcript_id\tgene_id\tn_exons\toverlapped_gene_count\toverlapped_gene_strand_type_count\toverlapped_gene_strands\toverlapped_genes"
    )?;

    let mut gtf_path = args.out.clone();
    gtf_path.add_extension("mark.gtf.gz");
    let gtf_file = File::create(&gtf_path)?;
    let mut gtf_writer = BufWriter::new(GzEncoder::new(gtf_file, Compression::default()));

    let mut reader = open_file_bufread(&args.input)?;
    let mut line = String::new();
    let mut line_no = 0usize;
    let mut tx_count = 0usize;
    let mut overlapped_tx_count = 0usize;

    while reader.read_line(&mut line)? != 0 {
        line_no += 1;
        let raw = line.trim_end();
        if raw.is_empty() || raw.starts_with('#') {
            writeln!(gtf_writer, "{raw}")?;
            line.clear();
            continue;
        }

        let fields: Vec<&str> = raw.split('\t').collect();
        if fields.len() < 9 {
            return Err(ToolError::FailedParseGTF {
                reason: format!("line {line_no}: expected 9 columns, got {}", fields.len()),
            });
        }
        if fields[2] != "transcript" {
            writeln!(gtf_writer, "{raw}")?;
            line.clear();
            continue;
        }

        tx_count += 1;
        let start = fields[3].parse::<u32>()?;
        let end = fields[4].parse::<u32>()?;
        let tx_id = parse_gtf_attr_value(fields[8], "transcript_id").ok_or_else(|| {
            ToolError::FailedParseGTF {
                reason: format!("line {line_no}: missing transcript_id"),
            }
        })?;
        let gene_id = parse_gtf_attr_value(fields[8], "gene_id").unwrap_or_default();
        let n_exons = parse_gtf_attr_value(fields[8], "ISOM_EXONS").unwrap_or("NA".to_string());

        let hits = gene_db.query_overlaps_range_all_strands(fields[0], start, end);
        let mut gene_hits = BTreeSet::new();
        let mut gene_ids = BTreeSet::new();
        let mut strands = BTreeSet::new();
        for hit in hits {
            strands.insert(strand_label(hit.strand));
            gene_ids.insert(hit.id.clone());
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

        write!(gtf_writer, "{}", fields[..8].join("\t"))?;
        writeln!(
            gtf_writer,
            "\t{}",
            append_gtf_attr(fields[8], "ISOM_OVLP_GENE", &format_ovlp_gene(&gene_ids)),
        )?;

        writeln!(
            tsv_writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            tx_id,
            gene_id,
            n_exons,
            gene_hits.len(),
            strands.len(),
            join_strs_or_na(strands),
            join_strings_or_na(gene_hits),
        )?;
        line.clear();
    }
    gtf_writer.flush()?;
    tsv_writer.flush()?;

    let stats = MarkStats {
        transcript_count: tx_count,
        overlapped_transcript_count: overlapped_tx_count,
    };
    print_json_block("Mark stats", &stats);

    let mut stats_path = args.out.clone();
    stats_path.add_extension("mark_stats.json");
    std::fs::write(&stats_path, serde_json::to_string_pretty(&stats)?)?;

    info!("GTF saved to: {}", gtf_path.display());
    info!("TSV saved to: {}", tsv_path.display());
    info!("Stats saved to: {}", stats_path.display());
    info!("Finished!");

    Ok(())
}

#[derive(Debug, Serialize)]
struct MarkStats {
    transcript_count: usize,
    overlapped_transcript_count: usize,
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

fn format_ovlp_gene(gene_ids: &BTreeSet<String>) -> String {
    if gene_ids.is_empty() {
        return "0:".to_string();
    }
    format!(
        "{}:{},",
        gene_ids.len(),
        gene_ids.iter().cloned().collect::<Vec<_>>().join(",")
    )
}

fn append_gtf_attr(attrs: &str, key: &str, value: &str) -> String {
    let attrs = attrs.trim_end().trim_end_matches(';').trim_end();
    format!("{attrs}; {key} \"{value}\";")
}

fn join_strings_or_na(values: BTreeSet<String>) -> String {
    if values.is_empty() {
        "NA".to_string()
    } else {
        values.into_iter().collect::<Vec<_>>().join(",")
    }
}
