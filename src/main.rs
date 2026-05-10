use std::path::PathBuf;

use clap::{ArgAction, Parser, Subcommand};
use log::error;
use serde::Serialize;
pub mod utils;
use utils::greetings2;
pub mod core;

use crate::{index::run_index, merge::run_merge};
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
        help = "Input transcript index files (.isomx)",
        required = true,
        num_args = 1..
    )]
    pub inputs: Vec<PathBuf>,

    #[clap(
        short = 'r',
        long = "refsites",
        help = "Reference splice sites in BED format for splice-site wobble comparison (optional)"
    )]
    pub refsites: Option<PathBuf>,

    #[clap(
        long = "reftss",
        help = "Reference TSS sites in BED format for splice-site wobble comparison (optional)"
    )]
    pub reftss: Option<PathBuf>,

    #[clap(
        long = "reftes",
        help = "Reference TES sites in BED format for splice-site wobble comparison (optional)"
    )]
    pub reftes: Option<PathBuf>,

    #[clap(
        short = 'D',
        long = "wob-d-nc",
        help = "Allowed donor-site wobble when attaching noncanonical transcripts to canonical backbones",
        default_value_t = 3
    )]
    pub wob_d_nc: u32,

    // Those
    #[clap(
        short = 'A',
        long = "wob-a-nc",
        help = "Allowed acceptor-site wobble when attaching noncanonical transcripts to canonical backbones",
        default_value_t = 3
    )]
    pub wob_a_nc: u32,

    #[clap(
        short = 'U',
        long = "wob-u-nc",
        visible_alias = "wob_u_nc",
        help = "Wobble for non canonical unstrand transcript",
        default_value_t = 3
    )]
    pub wob_u_nc: u32,

    #[clap(
        short = 'd',
        long = "wob-d",
        visible_alias = "wob_d",
        help = "Allowed donor-site wobble when merging canonical transcripts",
        default_value_t = 0
    )]
    pub wob_d: u32,

    #[clap(
        short = 'a',
        long = "wob-a",
        visible_alias = "wob_d",
        help = "Allowed acceptor-site wobble when merging canonical transcripts",
        default_value_t = 0
    )]
    pub wob_a: u32,

    #[clap(
        short = 'u',
        long = "wob-u",
        visible_alias = "wob_u",
        help = "Wobble for canonical unstrand transcript",
        default_value_t = 3
    )]
    pub wob_u: u32,

    #[clap(
        short = 's',
        long = "tss-wob",
        visible_alias = "tss_wob",
        help = "Wobble for merging tss in canonical transcirpts",
        default_value_t = 50
    )]
    pub tss_wob: u32,

    #[clap(
        short = 'e',
        long = "tes-wob",
        visible_alias = "tes_wob",
        help = "Wobble for merging tes in canonical transcirpts",
        default_value_t = 50
    )]
    pub tes_wob: u32,

    #[clap(
        short = 'S',
        long = "tss-wob-nc",
        visible_alias = "tss_wob_nc",
        help = "Wobble for merging tss in non canonical transcirpts",
        default_value_t = 50
    )]
    pub tss_wob_nc: u32,

    #[clap(
        short = 'E',
        long = "tes-wob-nc",
        visible_alias = "tes_wob_nc",
        help = "Wobble for merging tes in non canonical transcirpts",
        default_value_t = 50
    )]
    pub tes_wob_nc: u32,

    #[clap(
        short = 't',
        long = "terminal-merge",
        visible_alias = "terminal_merge",
        help = "Terminal merge mode for canonical transcirpts",
        value_enum,
        default_value = "both"
    )]
    pub terminal_merge: TerminalMergeMode,

    #[clap(
        short = 'T',
        long = "terminal-merge-nc",
        visible_alias = "terminal_merge_nc",
        help = "Terminal merge mode for non canonical transcirpts",
        value_enum,
        default_value = "both"
    )]
    pub terminal_merge_nc: TerminalMergeMode,

    #[clap(
        long = "tss-policy",
        help = "How to choose the representative TSS after transcripts have already been merged into one group",
        value_enum,
        default_value_t = MergePolicy::Major
    )]
    pub tss_policy: MergePolicy,

    #[clap(
        long = "tes-policy",
        help = "How to choose the representative TES after transcripts have already been merged into one group",
        value_enum,
        default_value_t = MergePolicy::Major
    )]
    pub tes_policy: MergePolicy,

    #[clap(
        long = "splice-policy",
        help = "How to choose the representative splice junction after transcripts have already been merged into one group",
        value_enum,
        default_value_t = MergePolicy::Major
    )]
    pub splice_policy: MergePolicy,

    #[clap(
        long = "mono-ovlp",
        help = "Minimum reciprocal overlap required for mono-exon merge",
        default_value_t = 0.9
    )]
    pub mono_ovlp: f64,

    #[clap(
        long = "mono-policy",
        help = "How to choose the representative start and end for monon exon transcripts have already been merged into one group",
        value_enum,
        default_value_t = MergePolicy::Major
    )]
    pub mono_policy: MergePolicy,

    #[clap(
        long = "unstranded-terminal-wobble",
        visible_alias = "uterminalwob",
        help = "Terminal wobble for unstranded transcripts",
        default_value_t = 50
    )]
    pub unstrand_terminal_wob: u32,

    #[clap(
        long = "sx-max",
        help = "Maximum exon length treated as a small-exon rescue candidate",
        default_value_t = 15
    )]
    pub sx_max: u32,

    #[clap(
        long = "junc-diff",
        help = "Maximum junction-count difference allowed for small-exon collapse rescue",
        default_value_t = 1
    )]
    pub junc_diff: u32,

    #[clap(
        long = "shift-rescue",
        help = "Enable rescue when only local small-exon boundaries shift",
        action = ArgAction::Set,
        default_value_t = true
    )]
    pub shift_rescue: bool,

    #[clap(
        long = "collapse-rescue",
        help = "Enable rescue when one or a few small exons collapse or disappear",
        action = ArgAction::Set,
        default_value_t = false
    )]
    pub collapse_rescue: bool,

    #[clap(
        long = "hash-rescue",
        help = "Require sequence or reference hash support before applying rescue rules",
        action = ArgAction::Set,
        default_value_t = true
    )]
    pub hash_rescue: bool,

    #[clap(short = 'o', long = "out", help = "Union output gtf file")]
    pub out: PathBuf,
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
