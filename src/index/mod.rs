use std::{
    fs::{self, File},
    io::Write,
    path::PathBuf,
};

use anyhow::{Context, bail};
use log::{error, info, warn};
use serde::Serialize;

use crate::{
    IndexArgs,
    traits::ArgValidate,
    utils::{greetings2, print_json_block, require_file, save_json_block},
};
pub use anyhow::Result as AnyResult;
use libgtf::{BuildConfig, BuildEvent, BuildReport, build_index_with_events};

#[derive(Debug, Serialize)]
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

impl From<BuildReport> for IndexStats {
    fn from(report: BuildReport) -> Self {
        Self {
            transcript_count: report.transcript_count,
            gene_count: report.gene_count,
            skipped_transcript_cnt: report.skipped_transcript_count,
            skipped_gene_cnt: report.skipped_gene_count,
            missing_seqid_cnt: report.missing_seqid_count,
            missing_seqids: report.missing_seqids,
            plus_strand_tx_cnt: report.plus_strand_transcript_count,
            minus_strand_tx_cnt: report.minus_strand_transcript_count,
            unknown_strand_tx_cnt: report.unknown_strand_transcript_count,
            mono_exon_tx_cnt: report.mono_exon_transcript_count,
            multi_exon_tx_cnt: report.multi_exon_transcript_count,
            all_canonical_tx_cnt: report.all_canonical_transcript_count,
            partial_canonical_tx_cnt: report.partial_canonical_transcript_count,
            non_canonical_tx_cnt: report.non_canonical_transcript_count,
            junction_cnt: report.junction_count,
            canonical_junction_cnt: report.canonical_junction_count,
            non_canonical_junction_cnt: report.non_canonical_junction_count,
            canonical_junction_ratio: report.canonical_junction_ratio,
        }
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

    let mut output_cleanup = OutputCleanup::new();

    let out_base = if let Some(out) = &args.out {
        out.clone()
    } else {
        args.input.clone()
    };

    let mut param_path = out_base.clone();
    param_path.add_extension("index_params.json");
    save_json_block(&param_path, &args)?;
    output_cleanup.track(param_path);

    if !args.quiet {
        info!("Creating isomatch index for {}", args.input.display());
    }

    let isomx_path = if let Some(out) = &args.out {
        let mut out_path = out.clone();
        out_path.add_extension("isomx");
        out_path
    } else {
        let mut default_out = args.input.clone();
        default_out.add_extension("isomx");
        default_out
    };

    let mut isoms_path = out_base.clone();
    isoms_path.add_extension("isoms");

    let config = BuildConfig {
        gtf_path: args.input.clone(),
        reference_fasta: args.ref_fa.clone(),
        transcript_fasta: args.seqfa.clone(),
        isomx_output: isomx_path.clone(),
        isoms_output: isoms_path.clone(),
        skip_missing_reference_seqids: args.skip_missing_ref_chr,
        temp_dir: None,
    };

    let quiet = args.quiet;
    let result = build_index_with_events(&config, |event| match event {
        BuildEvent::LoadingFasta if !quiet => {
            info!("Loading Reference and/or Sequence FASTA...");
        }
        BuildEvent::IndexingGtf if !quiet => info!("Indexing GTF"),
        BuildEvent::InitializingBuilders => {
            output_cleanup.track(isomx_path.clone());
            output_cleanup.track(isoms_path.clone());
            if !quiet {
                info!("Initializing Builder");
            }
        }
        BuildEvent::MissingReferenceSequence { seqid } => warn!(
            "Reference FASTA is missing seqid '{}'; transcripts on this seqid will be skipped",
            seqid
        ),
        BuildEvent::ProcessingChromosome { name } if !quiet => {
            info!("Processing chromosome {}", name);
        }
        BuildEvent::SkippingChromosome { name } if !quiet => {
            info!(
                "Skipping chromosome {} because it is absent from the reference FASTA",
                name
            );
        }
        BuildEvent::Finalizing => {}
        _ => {}
    });

    let report = match result {
        Ok(report) => report,
        Err(libgtf::error::Error::MissingReferenceSequences { seqids }) => {
            bail!(
                "Reference FASTA is missing {} seqid(s) required by the GTF: {}. Rerun with --skip-missing-ref-chr to warn and skip these transcripts. ",
                seqids.len(),
                seqids.join(", ")
            );
        }
        Err(libgtf::error::Error::NoIndexableSequences) => {
            bail!("No indexable seqids remain after filtering against the reference FASTA");
        }
        Err(err) => return Err(err).context("cannot build index"),
    };
    let stats = IndexStats::from(report);

    if !args.quiet {
        info!("Index isomx saved to {:?}", isomx_path);
        info!("Sidecar isoms saved to {:?}", isoms_path);
    }

    let mut index_info_path = out_base.clone();
    index_info_path.add_extension("index_info.json");
    let mut isomx_info_writer = File::create(&index_info_path)?;
    output_cleanup.track(index_info_path);
    if !args.quiet {
        print_json_block("Index summary", &stats);
    }

    let info_json = serde_json::to_string_pretty(&stats)?;

    isomx_info_writer.write_all(info_json.as_bytes())?;
    isomx_info_writer.flush()?;
    output_cleanup.disarm();

    if !args.quiet {
        info!("Finished!");
    }
    Ok(())
}
