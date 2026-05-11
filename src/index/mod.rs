use std::collections::HashSet;

use anyhow::{Context, bail};
use log::{error, info, warn};
use serde::Serialize;

use crate::{
    IndexArgs,
    core::tx_strand::ISOMSTRAND,
    fasta::{self, FastaReader},
    gtf::{self, profile_gtf},
    index::format::ChromBlockBuilder,
    traits::ArgValidate,
    utils::print_json_block,
};
pub use anyhow::Result;
pub mod builder;
pub mod format;
pub mod index_error;
pub mod reader;

#[derive(Debug, Default, Serialize)]
pub struct IndexStats {
    pub transcript_count: u64,
    pub gene_count: u64,
    pub skipped_transcript_count: u64,
    pub skipped_gene_count: u64,
    pub skipped_missing_ref_seqid_count: u64,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub skipped_missing_ref_seqids: Vec<String>,
    pub plus_strand_count: u64,
    pub minus_strand_count: u64,
    pub unknown_strand_count: u64,
    pub mono_exon_count: u64,
    pub multi_exon_count: u64,
    pub junction_count: u64,
    pub canonical_junction_count: u64,
    pub noncanonical_junction_count: u64,
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
            ISOMSTRAND::Minus => self.minus_strand_count += 1,
            ISOMSTRAND::Plus => self.plus_strand_count += 1,
            ISOMSTRAND::Unknown => self.unknown_strand_count += 1,
        }

        if exon_count <= 1 {
            self.mono_exon_count += 1;
            return;
        }

        self.multi_exon_count += 1;

        let junction_count = (exon_count - 1) as u64;
        let canonical_junction_count = canonical_junction_count as u64;

        self.junction_count += junction_count;
        self.canonical_junction_count += canonical_junction_count;
        self.noncanonical_junction_count += junction_count - canonical_junction_count;
    }

    pub fn observe_skipped_tx(&mut self, gene_id: &str) {
        self.skipped_transcript_count += 1;
        self.skipped_gene_ids.insert(gene_id.to_string());
    }

    pub fn note_skipped_ref_seqids(&mut self, seqids: Vec<String>) {
        self.skipped_missing_ref_seqid_count = seqids.len() as u64;
        self.skipped_missing_ref_seqids = seqids;
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
    let gtf_size = std::fs::metadata(&args.input)
        .expect("Can not read gtf file metadata")
        .len();
    let md5 = crate::utils::checksum_file(&args.input)
        .expect("Can not read gtf file for checksum")
        .0;
    let profiled_chrom_names = profile_gtf(&args.input)
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

    let output_path = if let Some(out) = &args.out {
        out.clone()
    } else {
        let mut default_out = args.input.clone();
        default_out.add_extension("isomx");
        default_out
    };

    info!("Initializing Builder...");
    let mut builder = builder::IndexBuilder::new(
        std::fs::File::create(&output_path).expect("Can not create output file"),
        chrom_names,
        gtf_size,
        md5,
        true,
        args.seqfa.is_some(),
    )
    .expect("Can not init index builder");

    info!("Indexing GTF...");
    let gtf_reader = gtf::MyGTFReader::new(&args.input)
        .with_context(|| format!("Can not open GTF file: {}", args.input.display()))?;

    let mut current_chrom = String::new();
    let mut chrom_id = 0u16;
    let mut chrom_block: Option<ChromBlockBuilder> = None;

    for tx_record in gtf_reader {
        if current_chrom != tx_record.chrom {
            if let Some(cb) = chrom_block.take() {
                builder.add_chrom(cb)?;
            }
            current_chrom = tx_record.chrom.clone();
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
        if missing_ref_seqid_set.contains(&tx_record.chrom) {
            stats.observe_skipped_tx(&tx_record.gene_id);
            continue;
        }
        chrom_block
            .as_mut()
            .expect("Can not access chromblock")
            .add_tx(tx_record, &mut ref_far, &mut seq_far, &mut stats)?;
    }

    if let Some(cb) = chrom_block.take() {
        builder.add_chrom(cb)?;
    }
    builder.finalize()?;
    stats.finalize();

    info!("Index written to {:?}", output_path);
    print_json_block("Index summary", &stats);
    info!("Fnished!");
    Ok(())
}
