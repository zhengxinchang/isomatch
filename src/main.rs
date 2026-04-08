use std::path::PathBuf;

use clap::{ArgAction, Parser, Subcommand, ValueEnum};
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
#[derive(Parser, Debug, Serialize, Clone)]
#[clap(
    name = "isomatch",
    version = "0.1.0",
    author = "xinchang zheng <zhengxc93@gmail.com>",
    about = "A verstile tool for isoform comparison and correction",
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
    Compare(CompareArgs),
    Annotate(AnnotateArgs),
}

#[derive(Copy, Clone, Debug, Serialize, ValueEnum)]
pub enum ReprEndPolicy {
    Outer,
    Inner,
}

#[derive(Parser, Debug, Serialize, Clone)]
#[clap(
    about = "index GTF files with sequence information, will be called by other subcommands, but can also be used independently
"
)]
pub struct IndexArgs {
    pub input: PathBuf,

    #[clap(short = 'r', long = "reffa", help = "Reference fasta file")]
    pub reffa: PathBuf,

    // hide this because the tx seqs usually not have paried sequence.
    #[clap(skip = None)]
    // #[clap(short = 's', long = "seqfa", help = "Transcirpt sequence file in fasta format")]
    pub seqfa: Option<PathBuf>,

    #[clap(
        short = 'o',
        long = "out",
        help = "index file, default is input file with .idx suffix"
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
        short = 'd',
        long = "dwob",
        help = "Allowed donor-site coordinate wobble during junction comparison",
        default_value_t = 3
    )]
    pub dwob: u32,

    #[clap(
        short = 'a',
        long = "awob",
        help = "Allowed acceptor-site coordinate wobble during junction comparison",
        default_value_t = 3
    )]
    pub awob: u32,

    #[clap(
        short = 's',
        long = "tss",
        help = "TSS tolerance for final merge decisions",
        default_value_t = 50
    )]
    pub tss: u32,

    #[clap(
        short = 'e',
        long = "tes",
        help = "TES tolerance for final merge decisions",
        default_value_t = 100
    )]
    pub tes: u32,

    #[clap(
        long = "mono-ovlp",
        help = "Minimum reciprocal overlap required for mono-exon merge",
        default_value_t = 0.9
    )]
    pub mono_ovlp: f64,

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

    #[clap(
        long = "repr-ends",
        help = "How representative transcript ends are chosen after merge",
        value_enum,
        default_value_t = ReprEndPolicy::Outer
    )]
    pub repr_ends: ReprEndPolicy,

    #[clap(short = 'o', long = "out", help = "Union output gtf file")]
    pub out: PathBuf,
}

#[derive(Parser, Debug, Serialize, Clone)]
#[clap(
    about = "compare multiple query transcript sets in GTF format and report summary statistics
"
)]
pub struct CompareArgs {
    #[clap(help = "GTF files to compare, at lest two files", required = true, num_args = 2..)]
    pub input: Vec<PathBuf>,

    #[clap(
        short = 'o',
        long = "out",
        help = "Output prefix for comparison results"
    )]
    pub out: PathBuf,
}

#[derive(Parser, Debug, Serialize, Clone)]
#[clap(about = "annotate query transcripts with reference annotation
")]
pub struct AnnotateArgs {
    pub input: PathBuf,

    #[clap(
        short = 'r',
        long = "annotation",
        help = "Reference annotation gtf file"
    )]
    pub annotation: PathBuf,

    #[clap(
        short = 'o',
        long = "out",
        help = "Output gtf file with classification results"
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
    let mut cli = Cli::parse();

    // set env logger level to info by default, can be overridden by RUST_LOG env var
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    match cli {
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
        Cli {
            command: Commands::Compare(args),
        } => {
            greetings2(&args);
        }
        Cli {
            command: Commands::Annotate(args),
        } => {
            greetings2(&args);
        }
    }
}
