use std::{collections::HashSet, fs::File, io::Write};

use anyhow::{Context, bail};
use log::{error, info, warn};
use serde::Serialize;

use crate::{
    IndexArgs,
    core::tx_strand::ISOMSTRAND,
    // fasta::{self, FastaReader},
    // gtf::{self, profile_gtf},
    index::format::ChromBlockBuilder,
    traits::ArgValidate,
    utils::{greetings2, print_json_block},
};
pub use anyhow::Result as AnyResult;
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

#[derive(Debug, Default, Serialize)]
pub struct IndexStats {
    pub transcript_count: u64,
    pub gene_count: u64,
    pub skipped_transcript_cnt: u64,
    pub skipped_gene_cnt: u64,
    pub missing_seqid_cnt: u64,
    pub missing_seqids: Vec<String>,
    pub plus_strand_tx_cnt: u64,
    pub minus_strand_tx_cnt: u64,
    pub unknown_strand_tx_cnt: u64,
    pub mono_exon_tx_cnt: u64,
    pub multi_exon_tx_cnt: u64,
    pub all_canonical_tx_cnt: u64,
    pub partial_canonical_tx_cnt: u64,
    pub non_canonical_tx_cnt: u64,
    pub junction_cnt: u64,
    pub canonical_junction_cnt: u64,
    pub non_canonical_junction_cnt: u64,
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
            ISOMSTRAND::Minus => self.minus_strand_tx_cnt += 1,
            ISOMSTRAND::Plus => self.plus_strand_tx_cnt += 1,
            ISOMSTRAND::Unknown => self.unknown_strand_tx_cnt += 1,
        }

        if exon_count <= 1 {
            self.mono_exon_tx_cnt += 1;
            return;
        }

        self.multi_exon_tx_cnt += 1;

        let junction_count = (exon_count - 1) as u64;
        let canonical_junction_count = canonical_junction_count as u64;

        if canonical_junction_count == junction_count {
            self.all_canonical_tx_cnt += 1;
        } else if canonical_junction_count == 0 {
            self.non_canonical_tx_cnt += 1;
        } else {
            self.partial_canonical_tx_cnt += 1;
        }

        self.junction_cnt += junction_count;
        self.canonical_junction_cnt += canonical_junction_count;
        self.non_canonical_junction_cnt += junction_count - canonical_junction_count;
    }

    pub fn observe_skipped_tx(&mut self, gene_id: &str) {
        self.skipped_transcript_cnt += 1;
        self.skipped_gene_ids.insert(gene_id.to_string());
    }

    pub fn note_skipped_ref_seqids(&mut self, seqids: Vec<String>) {
        self.missing_seqid_cnt = seqids.len() as u64;
        self.missing_seqids = seqids;
    }

    pub fn finalize(&mut self) {
        self.gene_count = self.gene_ids.len() as u64;
        self.skipped_gene_cnt = self.skipped_gene_ids.len() as u64;
        self.canonical_junction_ratio = if self.junction_cnt == 0 {
            0.0
        } else {
            self.canonical_junction_cnt as f64 / self.junction_cnt as f64
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

        if !self.ref_fa.exists() {
            error_msg.push_str(&format!(
                "\nReference FASTA file does not exist: {:?}",
                self.ref_fa
            ));
            has_error = true;
        }

        let mut fai1 = self.ref_fa.clone();
        fai1.add_extension("fai");
        if !fai1.exists() {
            error_msg.push_str(&format!(
                "\nReference FASTA index file does not exist: {:?}, use ' samtools faidx {} ' to create one.",
                fai1,
                self.ref_fa.display()
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

pub fn run_index(args: &mut IndexArgs) -> AnyResult<()> {
    if !args.quiet {
        greetings2(&args);
    }

    args.validate();
    let mut stats = IndexStats::default();

    if !args.quiet {
        info!("Creating isomatch index for {}", args.input.display());

        info!("Loading Reference and/or Sequence FASTA...");
    }

    let mut ref_far = FastaReader::open(args.ref_fa.clone(), fasta::FaType::Ref)
        .with_context(|| format!("Can not load reference sequence: {}", args.ref_fa.display()))?;

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

    if !args.quiet {
        info!("Profiling GTF");
    }
    let (profiled_chrom_names, md5, gtf_file_size) = profile_gtf(&args.input)?;

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
                "Reference FASTA is missing {} seqid(s) required by the GTF: {}. Rerun with --skip-missing-ref-chr to warn and skip these transcripts.",
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

    if !args.quiet {
        info!("Initializing Builder");
    }
    let missing_seqids_vec: Vec<String> = missing_ref_seqid_set.iter().cloned().collect();
    let mut builder = builder::IndexBuilder::new(
        std::fs::File::create(&isomx_path).expect("Can not create output file"),
        chrom_names,
        gtf_file_size,
        md5,
        true,
        args.seqfa.is_some(),
        missing_seqids_vec,
    )
    .expect("Can not init index builder");

    if !args.quiet {
        info!("Indexing GTF");
    }
    let mut gtf_reader = gtf::MyGTFReader::new(&args.input)
        .with_context(|| format!("Can not open GTF file: {}", args.input.display()))?;

    let mut current_chrom = String::new();
    let mut chrom_id = 0u16;
    let mut chrom_block: Option<ChromBlockBuilder> = None;
    let mut next_written_tx_idx = 0u32;
    loop {
        let Some(mut tx_structure) = gtf_reader.next()? else {
            break;
        };

        if current_chrom != tx_structure.chrom {
            if let Some(cb) = chrom_block.take() {
                builder.add_chrom(cb)?;
            }
            current_chrom = tx_structure.chrom.clone();
            if missing_ref_seqid_set.contains(&current_chrom) {
                if !args.quiet {
                    info!(
                        "Skipping chromosome {} because it is absent from the reference FASTA",
                        &current_chrom
                    );
                }
                chrom_block = None;
            } else {
                chrom_id += 1;
                chrom_block = Some(ChromBlockBuilder::init(chrom_id));
                if !args.quiet {
                    info!("Processing chromosome {}", &current_chrom);
                }
            }
        }
        if missing_ref_seqid_set.contains(&tx_structure.chrom) {
            stats.observe_skipped_tx(&tx_structure.gene_id);
            continue;
        }

        tx_structure.set_gidx(next_written_tx_idx);
        chrom_block
            .as_mut()
            .expect("Can not access chromblock")
            .add_tx(tx_structure, &mut ref_far, &mut seq_far, &mut stats)?;

        next_written_tx_idx += 1;
    }

    if let Some(cb) = chrom_block.take() {
        builder.add_chrom(cb)?;
    }
    // isom_src_cache_builder.finalize()?;
    builder.finalize()?;
    stats.finalize();
    if !args.quiet {
        info!("Index written to {:?}", isomx_path);
    }
    // second pass to build the sidecar file isoms

    {
        use std::io::BufRead;
        if !args.quiet {
            info!("Profiling attributes sidecar file");
        }
        let index_file = File::open(&isomx_path).with_context(|| {
            format!(
                "cannot reopen index for second pass: {}",
                isomx_path.display()
            )
        })?;
        let mut index_reader = reader::IndexReader::open(index_file, 0)
            .with_context(|| "cannot open IndexReader for second pass")?;
        let txid_index = index_reader
            .build_txid_index()
            .with_context(|| "cannot build txid index")?;

        let mut isoms_path = isomx_path.clone();
        isoms_path.set_extension("isoms");

        let mut attr_builder = attributes_index::AttrIndexBuilder::init(
            &isoms_path,
            next_written_tx_idx as usize,
            &md5,
        )
        .with_context(|| format!("cannot init AttrIndexBuilder at {}", isoms_path.display()))?;

        let gtf_lines = crate::utils::open_file_bufread(&args.input).with_context(|| {
            format!(
                "cannot reopen GTF for second pass: {}",
                args.input.display()
            )
        })?;
        for line in gtf_lines.lines() {
            let line = line.with_context(|| "error reading GTF line in second pass")?;
            if line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.splitn(9, '\t').collect();
            if fields.len() < 9 || fields[2] != "transcript" {
                continue;
            }
            let attr_str = fields[8];
            let Some(tx_id) = parse_gtf_attr_value(attr_str, "transcript_id") else {
                continue;
            };
            let Some(&tx_gidx) = txid_index.get(&tx_id) else {
                continue; // filtered out (e.g. missing ref seqid)
            };
            attr_builder
                .dump_attr(attr_str.as_bytes().to_vec(), tx_gidx)
                .with_context(|| format!("dump_attr failed for {}", tx_id))?;
        }

        attr_builder
            .finish()
            .with_context(|| format!("cannot finalize isoms at {}", isoms_path.display()))?;
        if !args.quiet {
            info!("Sidecar isoms written to {:?}", isoms_path);
        }
    }

    let mut isomx_info_path = isomx_path.clone();
    isomx_info_path.add_extension("info.json");
    let mut isomx_info_writer = File::create(isomx_info_path)?;
    if !args.quiet {
        print_json_block("Index summary", &stats);
    }

    let info_json = serde_json::to_string_pretty(&stats)?;

    isomx_info_writer.write(info_json.as_bytes())?;
    isomx_info_writer.flush()?;

    if !args.quiet {
        info!("Finished!");
    }
    Ok(())
}
