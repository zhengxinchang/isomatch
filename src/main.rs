use std::path::PathBuf;

use clap::{Parser, Subcommand};
use log::error;
use serde::Serialize;
pub mod utils;
use utils::greetings2;
pub mod core;

use crate::index::run_index;
pub mod fasta;
pub mod gtf;
pub mod index;
pub mod traits;
pub mod merge;
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
#[clap(about = "merge multiple query transcript sets into union GTF")]
pub struct MergeArgs {
    pub inputs: Vec<PathBuf>,

    #[clap(
        short = 'd',
        long = "donor-wobble-bp",
        help = "Wobble base pairs for 5' splice sites",
        default_value_t = 0
    )]
    pub donor_wobble_bp: u32,

    #[clap(
        short = 'a',
        long = "acceptor-wobble-bp",
        help = "Wobble base pairs for 3' splice sites",
        default_value_t = 0
    )]
    pub acceptor_wobble_bp: u32,

    #[clap(
        short = 's',
        long = "tss-tolerance",
        help = "Transcription start site tolerance",
        default_value_t = 0
    )]
    pub tss: u32,

    #[clap(
        short = 'e',
        long = "tes-tolerance",
        help = "Transcription end site tolerance",
        default_value_t = 0
    )]
    pub tes: u32,

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
