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
use crate::merge::policy::{MergePolicyArg, TerminalRefineMode};
#[derive(Parser, Debug, Serialize, Clone)]
#[clap(
    name = "isomatch",
    version = "0.1.0",
    author = "xinchang zheng <zhengxc93@gmail.com>",
    about = "Isoform comparison and correction tools",
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
    about = "Build an indexed transcript set from a GTF and reference FASTA.
"
)]
pub struct IndexArgs {
    #[clap(help = "Input GTF")]
    pub input: PathBuf,

    #[clap(short = 'r', long = "reffa", help = "Reference FASTA")]
    pub reffa: PathBuf,

    // hide this because the tx seqs usually not have paried sequence.
    #[clap(skip = None)]
    // #[clap(short = 's', long = "seqfa", help = "Transcirpt sequence file in fasta format")]
    pub seqfa: Option<PathBuf>,

    #[clap(
        short = 'o',
        long = "out",
        help = "Output index path; defaults to <input>.isomx"
    )]
    pub out: Option<PathBuf>,

    #[clap(
        long = "skip-missing-ref-chr",
        action = ArgAction::SetTrue,
        help = "Skip transcripts on seqids absent from the reference FASTA"
    )]
    pub skip_missing_ref_chr: bool,
}

#[derive(Parser, Debug, Serialize, Clone)]
#[clap(about = "Merge indexed transcript sets into a union GTF")]
pub struct MergeArgs {
    #[clap(
        help_heading = "Input",
        help = "Indexed transcript sets to merge",
        required = true,
        num_args = 1..
    )]
    pub inputs: Vec<PathBuf>,

    // #[clap(
    //     short='l',
    //     long="list",
    //     help_heading = "Input",
    //     help = "Input manifest file for multiple transcript sets to merge. override default input. Format path \t name(optional)",
    //     required = false,
    // )]

    // pub inputs_list: Option<PathBuf>,
    // #[clap(
    //     short = 'r',
    //     long = "refsites",
    //     help_heading = "Guide Merge",
    //     help = "Optional BED file of reference splice sites"
    // )]
    // pub refsites: Option<PathBuf>,
    #[clap(
        short = 'd',
        long = "wob-d",
        help_heading = "Canonical Transcript Merge",
        help = "Canonical donor wobble in bp",
        default_value_t = 0,
        value_name = "BP"
    )]
    pub wob_d: u32,

    #[clap(
        short = 'a',
        long = "wob-a",
        help_heading = "Canonical Transcript Merge",
        help = "Canonical acceptor wobble in bp",
        default_value_t = 0,
        value_name = "BP"
    )]
    pub wob_a: u32,

    #[clap(
        short = 'u',
        long = "wob-u",
        help_heading = "Canonical Transcript Merge",
        help = "Canonical unstranded splice wobble in bp",
        default_value_t = 3,
        value_name = "BP"
    )]
    pub wob_u: u32,

    #[clap(
        short = 't',
        long = "terminal-refine",
        help_heading = "Canonical Transcript Merge",
        help = "Canonical terminal refine mode: 
        which terminal sites (TSS/TES) to use as criteria 
        for distinguishing transcripts; works with tss-wob and tes-wob
        ",
        value_enum,
        default_value = "both"
    )]
    pub terminal_refine: TerminalRefineMode,

    #[clap(
        short = 's',
        long = "tss-wob",
        help_heading = "Canonical Transcript Merge",
        help = "Canonical TSS wobble in bp",
        default_value_t = 50,
        value_name = "BP"
    )]
    pub tss_wob: u32,

    #[clap(
        short = 'e',
        long = "tes-wob",
        help_heading = "Canonical Transcript Merge",
        help = "Canonical TES wobble in bp",
        default_value_t = 50,
        value_name = "BP"
    )]
    pub tes_wob: u32,

    #[clap(
        short = 'D',
        long = "wob-d-nc",
        help_heading = "Non-Canonical Transcript Merge",
        help = "Non-canonical donor wobble in bp",
        default_value_t = 3,
        value_name = "BP"
    )]
    pub wob_d_nc: u32,

    #[clap(
        short = 'A',
        long = "wob-a-nc",
        help_heading = "Non-Canonical Transcript Merge",
        help = "Non-canonical acceptor wobble in bp",
        default_value_t = 3,
        value_name = "BP"
    )]
    pub wob_a_nc: u32,

    #[clap(
        short = 'U',
        long = "wob-u-nc",
        help_heading = "Non-Canonical Transcript Merge",
        help = "Non-canonical unstranded splice wobble in bp",
        default_value_t = 3,
        value_name = "BP"
    )]
    pub wob_u_nc: u32,

    #[clap(
        short = 'T',
        long = "terminal-refine-nc",
        help_heading = "Non-Canonical Transcript Merge",
        help = "Non-canonical terminal refine mode: 
        which terminal sites (TSS/TES) to use as criteria 
        for distinguishing transcripts; works with tss-wob-nc and tes-wob-nc
        ",
        value_enum,
        default_value = "both"
    )]
    pub terminal_refine_nc: TerminalRefineMode,

    #[clap(
        short = 'S',
        long = "tss-wob-nc",
        help_heading = "Non-Canonical Transcript Merge",
        help = "Non-canonical TSS wobble in bp",
        default_value_t = 50,
        value_name = "BP"
    )]
    pub tss_wob_nc: u32,

    #[clap(
        short = 'E',
        long = "tes-wob-nc",
        help_heading = "Non-Canonical Transcript Merge",
        help = "Non-canonical TES wobble in bp",
        default_value_t = 50,
        value_name = "BP"
    )]
    pub tes_wob_nc: u32,

    #[clap(
        long = "mono-ovlp",
        help_heading = "Mono Exon Transcirpt Merge",
        help = "Minimum reciprocal overlap for mono-exon merge",
        default_value_t = 0.9,
        value_name = "FLOAT"
    )]
    pub mono_ovlp: f64,

    #[clap(
        long = "guide-tss",
        help_heading = "Representative Selection",
        help = "Guide TSS BED file"
    )]
    pub guide_tss: Option<PathBuf>,

    #[clap(
        long = "guide-tes",
        help_heading = "Representative Selection",
        help = "Guide TES BED file"
    )]
    pub guide_tes: Option<PathBuf>,

    #[clap(
        long = "guide-tss-flank",
        help_heading = "Representative Selection",
        help = "TSS evidence search flank in bp; used only with --guide-tss",
        default_value_t = 10,
        value_name = "BP"
    )]
    pub guide_tss_flank: u32,

    #[clap(
        long = "guide-tes-flank",
        help_heading = "Representative Selection",
        help = "TES evidence search flank in bp; used only with --guide-tes",
        default_value_t = 10,
        value_name = "BP"
    )]
    pub guide_tes_flank: u32,

    #[clap(
        long = "chrmap",
        help_heading = "Representative Selection",
        help = "Chromosome name map (UCSC -> Ensembl); use when inputs lack chr prefixes"
    )]
    pub chrmap: Option<PathBuf>,

    #[clap(
        long = "splice-policy",
        help_heading = "Representative Selection",
        help = "Representative splice-junction policy. 
    longer = longest exon span (thus shortest intron); 
    shorter = shortest exon span; 
    major = most frequent junction (falls back to longer on tie)
    ",
        value_enum,
        default_value_t = MergePolicyArg::Major
    )]
    pub splice_policy: MergePolicyArg,

    #[clap(
        long = "tss-policy",
        help_heading = "Representative Selection",
        help = "Representative TSS policy. 
    longer = most upstream TSS; 
    shorter = most downstream TSS; 
    major = most frequent TSS (falls back to union on tie)
    ",
        value_enum,
        default_value_t = MergePolicyArg::Major
    )]
    pub tss_policy: MergePolicyArg,

    #[clap(
        long = "tes-policy",
        help_heading = "Representative Selection",
        help = "Representative TES policy. 
    longer = most downstream TES; 
    shorter = most upstream TES; 
    major = most frequent TES (falls back to union on tie)
    ",
        value_enum,
        default_value_t = MergePolicyArg::Major
    )]
    pub tes_policy: MergePolicyArg,

    #[clap(
        long = "mono-policy",
        help_heading = "Representative Selection",
        help = "Representative mono-exon boundary policy. 
    union = widest span; 
    intersect = narrowest span; 
    major = most frequent span (falls back to union on tie)
    ",
        value_enum,
        default_value_t = MergePolicyArg::Major
    )]
    pub mono_policy: MergePolicyArg,

    #[clap(
        short = 'o',
        long = "out",
        help_heading = "Output",
        help = "Output union GTF"
    )]
    pub out: PathBuf,
    // #[clap(
    //     long = "sx-max",
    //     help_heading = "Other",
    //     help = "Max exon length for small-exon rescue",
    //     default_value_t = 15
    // )]
    // pub sx_max: u32,

    // #[clap(
    //     long = "junc-diff",
    //     help_heading = "Other",
    //     help = "Max junction-count difference for collapse rescue",
    //     default_value_t = 1
    // )]
    // pub junc_diff: u32,

    // #[clap(
    //     long = "shift-rescue",
    //     help_heading = "Other",
    //     help = "Enable local small-exon shift rescue",
    //     action = ArgAction::Set,
    //     default_value_t = true
    // )]
    // pub shift_rescue: bool,
}

#[derive(Parser, Debug, Serialize, Clone)]
#[clap(about = "Annotate query transcripts with a reference annotation
")]
pub struct AnnotateArgs {
    #[clap(help = "Query GTF")]
    pub input: PathBuf,

    #[clap(short = 'r', long = "annotation", help = "Reference annotation GTF")]
    pub annotation: PathBuf,

    #[clap(short = 'o', long = "out", help = "Output annotated GTF")]
    pub out: PathBuf,

    #[clap(
        short = 'c',
        long = "classification",
        help = "Classification mode: squant3, gffcompare, or both",
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
            if let Err(e) = run_index(&mut args) {
                error!("{}", e);
                std::process::exit(1);
            }
        }
        Cli {
            command: Commands::Merge(args),
        } => {
            greetings2(&args);
            if let Err(e) = run_merge(args) {
                error!("{}", e);
                std::process::exit(1);
            }
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
