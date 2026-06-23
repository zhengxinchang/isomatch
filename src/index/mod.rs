use std::{
    collections::{HashMap, HashSet},
    fs::{self, File},
    io::Write,
    path::PathBuf,
};

use anyhow::{Context, bail};
use log::{error, info, warn};
use serde::Serialize;

use crate::{
    IndexArgs,
    core::tx_strand::ISOMSTRAND,
    // fasta::{self, FastaReader},
    index::format::ChromBlockBuilder,
    traits::ArgValidate,
    utils::{greetings2, print_json_block, require_file},
};
pub use anyhow::Result as AnyResult;
use fasta::FastaReader;
pub mod attributes_index;
pub mod builder;
pub mod fasta;
pub mod format;
pub mod gtf;
pub mod index_error;
pub mod reader;

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

struct OutputCleanup {
    paths: Vec<PathBuf>,
    armed: bool,
}

impl OutputCleanup {
    fn new() -> Self {
        Self {
            paths: Vec::new(),
            armed: true,
        }
    }

    fn track(&mut self, path: PathBuf) {
        self.paths.push(path);
    }

    fn disarm(&mut self) {
        self.armed = false;
    }
}

impl Drop for OutputCleanup {
    fn drop(&mut self) {
        if self.armed {
            for path in &self.paths {
                let _ = fs::remove_file(path);
            }
        }
    }
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
        let mut error_msg = String::new();
        let mut has_error = false;

        require_file(
            "Input GTF file",
            &self.input,
            &mut error_msg,
            &mut has_error,
        );
        require_file(
            "Reference FASTA file",
            &self.ref_fa,
            &mut error_msg,
            &mut has_error,
        );

        let mut fai1 = self.ref_fa.clone();
        fai1.add_extension("fai");
        if !require_file(
            "Reference FASTA index file",
            &fai1,
            &mut error_msg,
            &mut has_error,
        ) && !fai1.exists()
        {
            error_msg.push_str(&format!(
                ", use ' samtools faidx {} ' to create one.",
                self.ref_fa.display()
            ));
        }

        if let Some(seqfa) = &self.seqfa {
            require_file("Sequence FASTA file", seqfa, &mut error_msg, &mut has_error);

            let mut seqfai1 = seqfa.clone();
            seqfai1.add_extension("fai");
            if !require_file(
                "Sequence FASTA index file",
                &seqfai1,
                &mut error_msg,
                &mut has_error,
            ) && !seqfai1.exists()
            {
                error_msg.push_str(&format!(
                    ", use ' samtools faidx {} ' to create one.",
                    seqfa.display()
                ));
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
        info!("Indexing GTF");
    }

    let mut gtf_reader = gtf::MyGTFReader::new(&args.input)
        .with_context(|| format!("Can not open GTF file: {}", args.input.display()))?;
    let profile = gtf_reader.profile().clone();

    let missing_ref_seqids: Vec<String> = profile
        .chrom_names
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
                "Reference FASTA is missing {} seqid(s) required by the GTF: {}. Rerun with --skip-missing-ref-chr to warn and skip these transcripts. ",
                missing_ref_seqids.len(),
                missing_ref_seqids.join(", ")
            );
        }
    }

    let missing_ref_seqid_set: HashSet<String> = missing_ref_seqids.into_iter().collect();
    let missing_ref_chrom_ids: HashSet<gtf::ChromID> = missing_ref_seqid_set
        .iter()
        .filter_map(|chrom| profile.chrom_name_to_id.get(chrom).copied())
        .collect();
    let chrom_names: Vec<String> = profile
        .chrom_names
        .iter()
        .filter(|chrom| !missing_ref_seqid_set.contains(*chrom))
        .cloned()
        .collect();
    let output_chrom_ids: HashMap<gtf::ChromID, u16> = profile
        .chrom_names
        .iter()
        .filter(|chrom| !missing_ref_seqid_set.contains(*chrom))
        .enumerate()
        .filter_map(|(out_idx, chrom)| {
            profile
                .chrom_name_to_id
                .get(chrom)
                .map(|profile_id| (*profile_id, (out_idx + 1) as u16))
        })
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
    let mut output_cleanup = OutputCleanup::new();
    let isomx_file = File::create(&isomx_path)
        .with_context(|| format!("Can not create output file: {}", isomx_path.display()))?;
    output_cleanup.track(isomx_path.clone());
    let mut builder = builder::IndexBuilder::new(
        isomx_file,
        chrom_names,
        profile.file_size,
        profile.md5,
        true,
        args.seqfa.is_some(),
        missing_seqids_vec,
    )
    .with_context(|| format!("Can not init index builder at {}", isomx_path.display()))?;

    let mut isoms_path = isomx_path.clone();
    isoms_path.set_extension("isoms");
    let total_indexable_tx = gtf_reader.transcript_count_excluding(&missing_ref_chrom_ids);
    let isoms_file = File::create(&isoms_path)
        .with_context(|| format!("cannot create isoms at {}", isoms_path.display()))?;
    output_cleanup.track(isoms_path.clone());
    let mut attr_builder =
        attributes_index::AttrIndexBuilder::new(isoms_file, total_indexable_tx, &profile.md5)
            .with_context(|| format!("cannot init AttrIndexBuilder at {}", isoms_path.display()))?;

    let mut current_chrom_id = 0u16;
    let mut chrom_block: Option<ChromBlockBuilder> = None;
    let mut next_written_tx_idx = 0u64;
    loop {
        let Some(mut tx_structure) = gtf_reader.next()? else {
            break;
        };

        let chrom_name = gtf_reader
            .chrom_name(tx_structure.chrom_id)
            .with_context(|| format!("invalid chrom_id {}", tx_structure.chrom_id))?
            .to_string();

        if current_chrom_id != tx_structure.chrom_id {
            if let Some(cb) = chrom_block.take() {
                builder.add_chrom(cb)?;
            }
            current_chrom_id = tx_structure.chrom_id;
            if missing_ref_chrom_ids.contains(&current_chrom_id) {
                if !args.quiet {
                    info!(
                        "Skipping chromosome {} because it is absent from the reference FASTA",
                        chrom_name
                    );
                }
                chrom_block = None;
            } else {
                let output_chrom_id = *output_chrom_ids
                    .get(&current_chrom_id)
                    .with_context(|| format!("missing output chrom_id for {chrom_name}"))?;
                chrom_block = Some(ChromBlockBuilder::init(output_chrom_id));
                if !args.quiet {
                    info!("Processing chromosome {}", chrom_name);
                }
            }
        }
        if missing_ref_chrom_ids.contains(&tx_structure.chrom_id) {
            stats.observe_skipped_tx(&tx_structure.gene_id);
            continue;
        }

        tx_structure.set_gidx(next_written_tx_idx);
        let attr_string = tx_structure.attr_string.clone();
        chrom_block
            .as_mut()
            .context("Can not access chromblock")?
            .add_tx(
                tx_structure,
                &chrom_name,
                &mut ref_far,
                &mut seq_far,
                &mut stats,
            )?;

        if let Some(attr_string) = attr_string {
            attr_builder
                .dump_attr(attr_string, next_written_tx_idx)
                .with_context(|| format!("dump_attr failed for tx_idx {}", next_written_tx_idx))?;
        }

        next_written_tx_idx = next_written_tx_idx
            .checked_add(1)
            .context("written transcript index exceeded u64")?;
    }

    if let Some(cb) = chrom_block.take() {
        builder.add_chrom(cb)?;
    }
    // isom_src_cache_builder.finalize()?;
    builder.finalize()?;
    stats.finalize();

    attr_builder
        .finish()
        .with_context(|| format!("cannot finalize isoms at {}", isoms_path.display()))?;

    if !args.quiet {
        info!("Index isomx saved to {:?}", isomx_path);
        info!("Sidecar isoms saved to {:?}", isoms_path);
    }

    let mut isomx_info_path = isomx_path.clone();
    isomx_info_path.add_extension("info.json");
    let mut isomx_info_writer = File::create(&isomx_info_path)?;
    output_cleanup.track(isomx_info_path);
    if !args.quiet {
        print_json_block("Index summary", &stats);
    }

    let info_json = serde_json::to_string_pretty(&stats)?;

    isomx_info_writer.write(info_json.as_bytes())?;
    isomx_info_writer.flush()?;
    output_cleanup.disarm();

    if !args.quiet {
        info!("Finished!");
    }
    Ok(())
}
