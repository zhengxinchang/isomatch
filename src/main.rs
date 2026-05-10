use std::path::PathBuf;

use clap::{ArgAction, Parser, Subcommand};
use log::error;
use serde::Serialize;
pub mod utils;
use utils::greetings2;
pub mod core;

use crate::{index::run_index, merge::run_merge};
pub mod constants;
pub mod fasta;
pub mod gtf;
pub mod index;
pub mod merge;
pub mod traits;
use crate::merge::policy::{MergePolicy, TerminalMergeMode};
#[derive(Parser, Debug, Serialize, Clone)]
#[clap(
    name = "isomatch",
    version = "0.1.0",
    author = "xinchang zheng <zhengxc93@gmail.com>",
    about = "A versatile tool for isoform comparison and correction",
    after_long_help = "

> [Examples]




"
)]
pub struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug, Serialize, Clone)]
pub enum Commands {
    Index(IndexArgs),
    Merge(MergeArgs),
    // Bench(BenchArgs),
    Annotate(AnnotateArgs),
}

#[derive(Parser, Debug, Serialize, Clone)]
#[clap(
    about = "Index GTF files with sequence information. This command is used by other subcommands, but can also be run independently.
"
)]
pub struct IndexArgs {
    #[clap(help = "Input GTF file")]
    pub input: PathBuf,

    #[clap(short = 'r', long = "reffa", help = "Reference FASTA file")]
    pub reffa: PathBuf,

    // hide this because the tx seqs usually not have paried sequence.
    #[clap(skip = None)]
    // #[clap(short = 's', long = "seqfa", help = "Transcirpt sequence file in fasta format")]
    pub seqfa: Option<PathBuf>,

    #[clap(
        short = 'o',
        long = "out",
        help = "Output index file path; defaults to the input path with an .isomx suffix"
    )]
    pub out: Option<PathBuf>,
}

#[derive(Parser, Debug, Serialize, Clone)]
#[clap(about = "merge multiple indexed transcript sets into a union GTF")]
pub struct MergeArgs {
    #[clap(
        help_heading = "Input",
        help = "Input transcript sets to merge",
        required = true,
        num_args = 1..
    )]
    pub inputs: Vec<PathBuf>,

    #[clap(
        short = 'r',
        long = "refsites",
        help_heading = "Guide Merge",
        help = "Optional BED file of reference splice sites"
    )]
    pub refsites: Option<PathBuf>,

    #[clap(
        long = "reftss",
        help_heading = "Guide Merge",
        help = "Optional BED file of reference TSS sites"
    )]
    pub reftss: Option<PathBuf>,

    #[clap(
        long = "reftes",
        help_heading = "Guide Merge",
        help = "Optional BED file of reference TES sites"
    )]
    pub reftes: Option<PathBuf>,

    #[clap(
        short = 'd',
        long = "wob-d",
        help_heading = "Canonical Transcript Merge",
        help = "Donor wobble for canonical splice-junction merge",
        default_value_t = 0
    )]
    pub wob_d: u32,

    #[clap(
        short = 'a',
        long = "wob-a",
        help_heading = "Canonical Transcript Merge",
        help = "Acceptor wobble for canonical splice-junction merge",
        default_value_t = 0
    )]
    pub wob_a: u32,

    #[clap(
        short = 'u',
        long = "wob-u",
        help_heading = "Canonical Transcript Merge",
        help = "Wobble for unstranded canonical splice-junction merge",
        default_value_t = 3
    )]
    pub wob_u: u32,

    #[clap(
        short = 's',
        long = "tss-wob",
        help_heading = "Canonical Transcript Merge",
        help = "TSS wobble for canonical terminal refinement",
        default_value_t = 50
    )]
    pub tss_wob: u32,

    #[clap(
        short = 'e',
        long = "tes-wob",
        help_heading = "Canonical Transcript Merge",
        help = "TES wobble for canonical terminal refinement",
        default_value_t = 50
    )]
    pub tes_wob: u32,

    #[clap(
        short = 't',
        long = "terminal-merge",
        help_heading = "Canonical Transcript Merge",
        help = "How canonical groups are refined by TSS/TES",
        value_enum,
        default_value = "both"
    )]
    pub terminal_merge: TerminalMergeMode,

    #[clap(
        short = 'D',
        long = "wob-d-nc",
        help_heading = "Non-Canonical Transcript Merge",
        help = "Donor wobble for non-canonical attachment/merge",
        default_value_t = 3
    )]
    pub wob_d_nc: u32,

    #[clap(
        short = 'A',
        long = "wob-a-nc",
        help_heading = "Non-Canonical Transcript Merge",
        help = "Acceptor wobble for non-canonical attachment/merge",
        default_value_t = 3
    )]
    pub wob_a_nc: u32,

    #[clap(
        short = 'U',
        long = "wob-u-nc",
        help_heading = "Non-Canonical Transcript Merge",
        help = "Wobble for unstranded non-canonical splice-junction merge",
        default_value_t = 3
    )]
    pub wob_u_nc: u32,

    #[clap(
        short = 'S',
        long = "tss-wob-nc",
        help_heading = "Non-Canonical Transcript Merge",
        help = "TSS wobble for non-canonical terminal refinement",
        default_value_t = 50
    )]
    pub tss_wob_nc: u32,

    #[clap(
        short = 'E',
        long = "tes-wob-nc",
        help_heading = "Non-Canonical Transcript Merge",
        help = "TES wobble for non-canonical terminal refinement",
        default_value_t = 50
    )]
    pub tes_wob_nc: u32,

    #[clap(
        short = 'T',
        long = "terminal-merge-nc",
        help_heading = "Non-Canonical Transcript Merge",
        help = "How non-canonical groups are refined by TSS/TES",
        value_enum,
        default_value = "both"
    )]
    pub terminal_merge_nc: TerminalMergeMode,

    #[clap(
        long = "splice-policy",
        help_heading = "Representative Selection",
        help = "How to choose the representative splice junction",
        value_enum,
        default_value_t = MergePolicy::Major
    )]
    pub splice_policy: MergePolicy,

    #[clap(
        long = "tss-policy",
        help_heading = "Representative Selection",
        help = "How to choose the representative TSS",
        value_enum,
        default_value_t = MergePolicy::Major
    )]
    pub tss_policy: MergePolicy,

    #[clap(
        long = "tes-policy",
        help_heading = "Representative Selection",
        help = "How to choose the representative TES",
        value_enum,
        default_value_t = MergePolicy::Major
    )]
    pub tes_policy: MergePolicy,

    #[clap(
        long = "mono-policy",
        help_heading = "Representative Selection",
        help = "How to choose the representative mono-exon boundary pair",
        value_enum,
        default_value_t = MergePolicy::Major
    )]
    pub mono_policy: MergePolicy,

    #[clap(
        short = 'o',
        long = "out",
        help_heading = "Output",
        help = "Output union GTF path"
    )]
    pub out: PathBuf,

    #[clap(
        long = "mono-ovlp",
        help_heading = "Other",
        help = "Minimum reciprocal overlap for mono-exon merge",
        default_value_t = 0.9
    )]
    pub mono_ovlp: f64,

    #[clap(
        long = "sx-max",
        help_heading = "Other",
        help = "Maximum exon length considered a small-exon rescue target",
        default_value_t = 15
    )]
    pub sx_max: u32,

    #[clap(
        long = "junc-diff",
        help_heading = "Other",
        help = "Maximum junction-count difference for collapse rescue",
        default_value_t = 1
    )]
    pub junc_diff: u32,

    #[clap(
        long = "shift-rescue",
        help_heading = "Other",
        help = "Enable rescue for local small-exon boundary shifts",
        action = ArgAction::Set,
        default_value_t = true
    )]
    pub shift_rescue: bool,
}

#[derive(Parser, Debug, Serialize, Clone)]
#[clap(about = "Annotate query transcripts with a reference annotation
")]
pub struct AnnotateArgs {
    #[clap(help = "Input GTF file to annotate")]
    pub input: PathBuf,

    #[clap(
        short = 'r',
        long = "annotation",
        help = "Reference annotation GTF file"
    )]
    pub annotation: PathBuf,

    #[clap(
        short = 'o',
        long = "out",
        help = "Output GTF file with classification results"
    )]
    pub out: PathBuf,

    #[clap(
        short = 'c',
        long = "classification",
        help = "Classification system to use (squant3, gffcompare, both), default is 'both'",
        default_value = "both"
    )]
    pub classification: String,
}

fn main() {
    // set env logger level to info by default, can be overridden by RUST_LOG env var
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    match Cli::parse() {
        Cli {
            command: Commands::Index(mut args),
        } => {
            greetings2(&args);
            run_index(&mut args)
                .map_err(|e| {
                    error!("{}", e);
                })
                .expect("Exit...");
        }
        Cli {
            command: Commands::Merge(args),
        } => {
            greetings2(&args);
            run_merge(args)
                .map_err(|e| {
                    error!("{}", e);
                })
                .expect("Exit...");
        }
        // Cli {
        //     command: Commands::Bench(args),
        // } => {
        //     greetings2(&args);
        // }
        Cli {
            command: Commands::Annotate(args),
        } => {
            greetings2(&args);
        }
    }
}
