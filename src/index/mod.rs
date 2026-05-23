use std::{collections::HashSet, fs::File, io::Write};

use anyhow::{Context, bail};
use log::{error, info, warn};
use serde::Serialize;

use crate::{
    IndexArgs,
    core::tx_strand::ISOMSTRAND,
    // fasta::{self, FastaReader},
    // gtf::{self, profile_gtf},
    index::{attributes_index::IsomAttrCacheBuilder, format::ChromBlockBuilder},
    traits::ArgValidate,
    utils::print_json_block,
};
pub use anyhow::Result;
use fasta::FastaReader;
use gtf::profile_gtf;
pub mod attributes_index;
pub mod builder;
pub mod fasta;
pub mod format;
pub mod gtf;
pub mod index_error;
pub mod reader;

fn parse_gtf_attr_value(attrs: &str, key: &str) -> Option<String> {
    attrs.split(';').find_map(|attr| {
        let attr = attr.trim();
        if !attr.starts_with(key) {
            return None;
        }

        if let Some(q_start) = attr.find('"') {
            let rest = &attr[q_start + 1..];
            let q_len = rest.find('"')?;
            return Some(rest[..q_len].to_string());
        }

        attr.split_ascii_whitespace()
            .nth(1)
            .map(ToString::to_string)
    })
}

fn parse_isom_tx_id(tx_id: &str) -> Option<u32> {
    tx_id
        .strip_prefix("ISOMT_")
        .and_then(|s| s.parse::<u32>().ok())
        .and_then(|id| id.checked_sub(1))
}

fn parse_isom_src_attr(tx_attr: &gtf::TxAttrs) -> Result<Option<(u32, String)>> {
    let attrs = tx_attr.attr_string();
    let Some(isom_src_vec) = parse_gtf_attr_value(attrs, "ISOM_SRC") else {
        return Ok(None);
    };

    let tx_id_str = parse_gtf_attr_value(attrs, "transcript_id")
        .context("merged transcript with ISOM_SRC is missing transcript_id")?;
    let tx_id = parse_isom_tx_id(&tx_id_str).with_context(|| {
        format!(
            "merged transcript_id '{}' does not follow the expected ISOMT_<n> format",
            tx_id_str
        )
    })?;

    Ok(Some((tx_id, isom_src_vec)))
}

#[derive(Debug, Default, Serialize)]
pub struct IndexStats {
    pub transcript_count: u64,
    pub gene_count: u64,
    pub skipped_transcript_count: u64,
    pub skipped_gene_count: u64,
    pub missing_seqid_count: u64,
    // #[serde(skip_serializing_if = "Vec::is_empty")]
    pub missing_seqids: Vec<String>,
    pub plus_strand_tx_count: u64,
    pub minus_strand_tx_count: u64,
    pub unknown_strand_tx_count: u64,
    pub mono_exon_tx_count: u64,
    pub multi_exon_tx_count: u64,
    pub all_canonical_tx_count: u64,
    pub partial_canonical_tx_count: u64,
    pub non_canonical_tx_count: u64,
    pub junction_count: u64,
    pub canonical_junction_count: u64,
    pub non_canonical_junction_count: u64,
    pub canonical_junction_ratio: f64,
    #[serde(skip_serializing)]
    gene_ids: HashSet<String>,
    #[serde(skip_serializing)]
    skipped_gene_ids: HashSet<String>,
}

impl IndexStats {
    pub fn observe_tx(
        &mut self,
        strand: ISOMSTRAND,
        exon_count: usize,
        canonical_junction_count: usize,
        gene_id: &str,
    ) {
        self.transcript_count += 1;
        self.gene_ids.insert(gene_id.to_string());

        match strand {
            ISOMSTRAND::Minus => self.minus_strand_tx_count += 1,
            ISOMSTRAND::Plus => self.plus_strand_tx_count += 1,
            ISOMSTRAND::Unknown => self.unknown_strand_tx_count += 1,
        }

        if exon_count <= 1 {
            self.mono_exon_tx_count += 1;
            return;
        }

        self.multi_exon_tx_count += 1;

        let junction_count = (exon_count - 1) as u64;
        let canonical_junction_count = canonical_junction_count as u64;

        if canonical_junction_count == junction_count {
            self.all_canonical_tx_count += 1;
        } else if canonical_junction_count == 0 {
            self.non_canonical_tx_count += 1;
        } else {
            self.partial_canonical_tx_count += 1;
        }

        self.junction_count += junction_count;
        self.canonical_junction_count += canonical_junction_count;
        self.non_canonical_junction_count += junction_count - canonical_junction_count;
    }

    pub fn observe_skipped_tx(&mut self, gene_id: &str) {
        self.skipped_transcript_count += 1;
        self.skipped_gene_ids.insert(gene_id.to_string());
    }

    pub fn note_skipped_ref_seqids(&mut self, seqids: Vec<String>) {
        self.missing_seqid_count = seqids.len() as u64;
        self.missing_seqids = seqids;
    }

    pub fn finalize(&mut self) {
        self.gene_count = self.gene_ids.len() as u64;
        self.skipped_gene_count = self.skipped_gene_ids.len() as u64;
        self.canonical_junction_ratio = if self.junction_count == 0 {
            0.0
        } else {
            self.canonical_junction_count as f64 / self.junction_count as f64
        };
    }
}

impl ArgValidate for IndexArgs {
    fn validate(&self) {
        let mut error_msg = "".to_string();
        let mut has_error = false;

        if !self.input.exists() {
            error_msg.push_str(&format!(
                "\nInput GTF file does not exist: {:?}",
                self.input
            ));
            has_error = true;
        }

        if !self.reffa.exists() {
            error_msg.push_str(&format!(
                "\nReference FASTA file does not exist: {:?}",
                self.reffa
            ));
            has_error = true;
        }

        let mut fai1 = self.reffa.clone();
        fai1.add_extension("fai");
        if !fai1.exists() {
            error_msg.push_str(&format!(
                "\nReference FASTA index file does not exist: {:?}, use ' samtools faidx {} ' to create one.",
                fai1,
                self.reffa.display()
            ));
            has_error = true;
        }

        if let Some(seqfa) = &self.seqfa {
            if !seqfa.exists() {
                error_msg.push_str(&format!(
                    "\nSequence FASTA file does not exist: {:?}",
                    seqfa
                ));
                has_error = true;
            }

            let mut seqfai1 = seqfa.clone();
            seqfai1.add_extension("fai");
            if !seqfai1.exists() {
                error_msg.push_str(&format!(
                    "\nSequence FASTA index for {:?} does not exist, use ' samtools faidx {} ' to create one.",
                    seqfai1,
                    seqfa.display()
                ));
                has_error = true;
            }
        }

        if has_error {
            error!("Error validating arguments: {}", error_msg);
            std::process::exit(1);
        }
    }
}

pub fn run_index(args: &mut IndexArgs) -> Result<()> {
    args.validate();
    let mut stats = IndexStats::default();

    info!("Creating isomatch index for {}", args.input.display());

    info!("Loading Reference and/or Sequence FASTA...");

    let mut ref_far = FastaReader::open(args.reffa.clone(), fasta::FaType::Ref)
        .with_context(|| format!("Can not load reference sequence: {}", args.reffa.display()))?;

    let mut seq_far = if let Some(seqfa) = &args.seqfa {
        Some(
            FastaReader::open(seqfa.clone(), fasta::FaType::Seq).with_context(|| {
                format!(
                    "Can not load sequence from reference genome: {}",
                    seqfa.display()
                )
            })?,
        )
    } else {
        None
    };

    info!("Profiling GTF...");
    let (profiled_chrom_names, md5, gtf_size) = profile_gtf(&args.input)
        .with_context(|| format!("Can not profile GTF file: {}", args.input.display()))?;

    let missing_ref_seqids: Vec<String> = profiled_chrom_names
        .iter()
        .filter(|chrom| !ref_far.contains(chrom))
        .cloned()
        .collect();

    if !missing_ref_seqids.is_empty() {
        if args.skip_missing_ref_chr {
            for seqid in &missing_ref_seqids {
                warn!(
                    "Reference FASTA is missing seqid '{}'; transcripts on this seqid will be skipped",
                    seqid
                );
            }
            stats.note_skipped_ref_seqids(missing_ref_seqids.clone());
        } else {
            bail!(
                "Reference FASTA is missing {} seqid(s) required by the GTF: {}. Rerun with --skip-missing-ref-seqids to warn and skip these transcripts.",
                missing_ref_seqids.len(),
                missing_ref_seqids.join(", ")
            );
        }
    }

    let missing_ref_seqid_set: HashSet<String> = missing_ref_seqids.into_iter().collect();
    let chrom_names: Vec<String> = profiled_chrom_names
        .into_iter()
        .filter(|chrom| !missing_ref_seqid_set.contains(chrom))
        .collect();

    if chrom_names.is_empty() {
        bail!("No indexable seqids remain after filtering against the reference FASTA");
    }

    let isomx_path = if let Some(out) = &args.out {
        out.clone()
    } else {
        let mut default_out = args.input.clone();
        default_out.add_extension("isomx");
        default_out
    };

    info!("Initializing Builder...");
    let mut builder = builder::IndexBuilder::new(
        std::fs::File::create(&isomx_path).expect("Can not create output file"),
        chrom_names,
        gtf_size,
        md5,
        true,
        args.seqfa.is_some(),
    )
    .expect("Can not init index builder");

    info!("Indexing GTF...");
    let mut gtf_reader = gtf::MyGTFReader::new(&args.input)
        .with_context(|| format!("Can not open GTF file: {}", args.input.display()))?;

    let mut current_chrom = String::new();
    let mut chrom_id = 0u16;
    let mut chrom_block: Option<ChromBlockBuilder> = None;
    let mut next_written_tx_idx: u32 = 0;

    // init isomsrccache

    let mut isom_src_cache = IsomAttrCacheBuilder::init(&isomx_path);

    // for mut tx_record in gtf_reader {
    loop {
        let Some(gtf_record) = gtf_reader.next()? else {
            break;
        };

        match gtf_record {
            gtf::GTFRecord::TxAttrs(tx_attr) => {
                if let Some((tx_id, isom_src_vec)) = parse_isom_src_attr(&tx_attr)? {
                    // We first store the span by the merged transcript id parsed
                    // from `transcript_id`. After structure indexing decides
                    // whether this transcript survives chromosome filtering, we
                    // remap that span onto the dense `TxBase.tx_idx`.
                    isom_src_cache.dump_isom_src_string(isom_src_vec, tx_id)?;
                }
            }
            gtf::GTFRecord::TxStructure(mut tx_structure) => {
                // process the tx_strucuture as previous
                if current_chrom != tx_structure.chrom {
                    if let Some(cb) = chrom_block.take() {
                        builder.add_chrom(cb)?;
                    }
                    current_chrom = tx_structure.chrom.clone();
                    if missing_ref_seqid_set.contains(&current_chrom) {
                        info!(
                            "Skipping chromosome {} because it is absent from the reference FASTA",
                            &current_chrom
                        );
                        chrom_block = None;
                    } else {
                        chrom_id += 1;
                        chrom_block = Some(ChromBlockBuilder::init(chrom_id));
                        info!("Processing chromosome {}", &current_chrom);
                    }
                }
                if missing_ref_seqid_set.contains(&tx_structure.chrom) {
                    stats.observe_skipped_tx(&tx_structure.gene_id);
                    continue;
                }
                // resign the no skip continue id for tx
                // just incase the missing ref seqid skip
                // the tx_idx in TxBase will be used to prject the txbase_idx to src records
                // in isoms file.
                if let Some(original_tx_id) = parse_isom_tx_id(&tx_structure.tx_id) {
                    isom_src_cache.project_tx_id(original_tx_id, next_written_tx_idx);
                }
                tx_structure.set_idx(next_written_tx_idx);
                chrom_block
                    .as_mut()
                    .expect("Can not access chromblock")
                    .add_tx(tx_structure, &mut ref_far, &mut seq_far, &mut stats)?;
                next_written_tx_idx += 1;
            }
        }
    }

    if let Some(cb) = chrom_block.take() {
        builder.add_chrom(cb)?;
    }
    isom_src_cache.finalize()?;
    builder.finalize()?;
    stats.finalize();

    info!("Index written to {:?}", isomx_path);

    let mut isomx_info_path = isomx_path.clone();
    isomx_info_path.add_extension("index_info.json");
    let mut isomx_info_writer = File::create(isomx_info_path)?;
    print_json_block("Index summary", &stats);

    let info_json = serde_json::to_string_pretty(&stats)?;

    isomx_info_writer.write(info_json.as_bytes())?;
    isomx_info_writer.flush()?;

    info!("Fnished!");
    Ok(())
}
