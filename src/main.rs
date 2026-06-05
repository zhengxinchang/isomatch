use std::path::PathBuf;

use clap::{ArgAction, Parser, Subcommand};
use log::error;
use serde::Serialize;
pub mod core;
pub mod utils;
use crate::classify::run_classify;
use crate::tools::chop::run_chop;
use crate::{index::run_index, merge::run_merge};
use clap::ValueEnum;
pub mod classify;
pub mod constants;
// pub mod fasta;
// pub mod gtf;
pub mod index;
pub mod merge;
pub mod tools;
pub mod traits;
use crate::merge::policy::{MergePolicyArg, TerminalRefineMode};
#[derive(Parser, Debug, Serialize, Clone)]
#[clap(
    name = "isomatch",
    version = env!("CARGO_PKG_VERSION"),
    author = "xinchang zheng <zhengxc93@gmail.com>",
    about = "A tool for improved transcript merging and classification",
    after_long_help = "

Author: Xinchang Zheng <zhengxc93@gmail.com>

Repository: https://github.com/zhengxinchang/isomatch

> [Examples]

Build indexes:

isomatch index --ref-fa ref.fa sample1.gtf.gz
isomatch index --ref-fa ref.fa sample2.gtf.gz

Merge transcript sets:

isomatch merge --ref-fa ref.fa -o merged sample1.gtf.gz sample2.gtf.gz sample3.gtf.gz
# outputs: merged.merged.gtf.gz  merged.track.tsv.gz  merged.merged_info.json

Merge with guide-based terminal selection (human grch38 can be found at repo):

isomatch merge --ref-fa ref.fa -o merged --guide-tss human.grch38.tss.bed --guide-tes human.grch38.tes.bed sample1.gtf.gz sample2.gtf.gz sample3.gtf.gz

Merge with wobble splice junction matching:

isomatch merge --ref-fa ref.fa -o merged -d 3 -a 3 -u 3 -D 5 -A 5 -U 5 sample1.gtf.gz sample2.gtf.gz sample3.gtf.gz

Classify query transcripts against a reference annotation:

isomatch classify --ref-fa ref.fa --ref-gtf reference.gtf.gz -o query_vs_ref query.gtf.gz
# outputs: query_vs_ref.classification.txt.gz  query_vs_ref.annotated.gtf.gz  query_vs_ref.classify_info.json

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
    Classify(ClassifyArgs),
    #[command(subcommand)]
    Tools(ToolssArgs),
}

#[derive(Parser, Debug, Serialize, Clone)]
#[clap(
    about = "Build an indexed transcript set from a GTF and matched reference FASTA.
"
)]
pub struct IndexArgs {
    #[clap(help_heading = "Input", help = "Input GTF")]
    pub input: PathBuf,

    #[clap(
        short = 'r',
        long = "ref-fa",
        help_heading = "Input",
        help = "Reference FASTA"
    )]
    pub ref_fa: PathBuf,

    // hide this because the tx seqs usually not have paried sequence.
    #[clap(skip = None)]
    // #[clap(short = 's', long = "seqfa", help = "Transcirpt sequence file in fasta format")]
    pub seqfa: Option<PathBuf>,

    #[clap(
        short = 'o',
        long = "out",
        help_heading = "Output",
        help = "Output indexes path; defaults to <input>.isomx and <input>.isoms"
    )]
    pub out: Option<PathBuf>,

    #[clap(
        long = "skip-missing-ref-chr",
        action = ArgAction::SetTrue,
        help_heading = "Other",
        help = "Skip transcripts on seqids absent from the reference FASTA"
    )]
    pub skip_missing_ref_chr: bool,

    #[clap(
        short= 'q',
        long = "quiet",
        help_heading = "Other",
        action = ArgAction::SetTrue,
        help = "Suppress non-warning messages; only warnings and errors are shown"
    )]
    pub quiet: bool,
}

#[derive(Parser, Debug, Serialize, Clone)]
#[clap(about = "Merge indexed transcript sets into a union GTF
")]
pub struct MergeArgs {
    #[clap(
        help_heading = "Input",
        help = "Indexed transcript sets to merge",
        required = true,
        num_args = 1..
    )]
    pub inputs: Vec<PathBuf>,

    #[clap(
        short = 'o',
        long = "out",
        help_heading = "Output",
        help = "Output prefix "
    )]
    pub out: PathBuf,

    #[clap(
        short = 'r',
        long = "ref-fa",
        help_heading = "Auto Indexing",
        help = "Reference FASTA"
    )]
    pub ref_fa: PathBuf,

    #[clap(
        long = "skip-missing-ref-chr",
        action = ArgAction::SetTrue,
        help_heading = "Auto Indexing",
        help = "Skip transcripts on seqids absent from the reference FASTA"
    )]
    pub skip_missing_ref_chr: bool,

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
}

#[derive(Parser, Debug, Serialize, Clone)]
#[clap(about = "Classify query transcripts with a reference annotation GTF
")]
pub struct ClassifyArgs {
    #[clap(help_heading = "Input", help = "Query GTF")]
    pub input: PathBuf,

    #[clap(
        short = 'o',
        long = "out",
        help_heading = "Output",
        help = "Output prefix"
    )]
    pub out: PathBuf,

    #[clap(
        short = 'g',
        long = "ref-gtf",
        help_heading = "Reference data",
        help = "Reference annotation GTF"
    )]
    pub ref_gtf: PathBuf,

    #[clap(
        short = 'r',
        long = "ref-fa",
        help_heading = "Reference data",
        help = "Reference FASTA"
    )]
    pub ref_fa: PathBuf,

    #[clap(
        short = 's',
        long = "guide-tss",
        help_heading = "Reference data",
        help = "Reference TSS regions in bed format"
    )]
    pub guide_tss: Option<PathBuf>,

    #[clap(
        long = "chrmap",
        help_heading = "Reference data",
        help = "Chromosome name map (UCSC -> Ensembl); use when inputs lack chr prefixes"
    )]
    pub chrmap: Option<PathBuf>,

    #[clap(
        short = 'e',
        long = "guide-tes",
        help_heading = "Reference data",
        help = "Reference TES regions in bed format"
    )]
    pub guide_tes: Option<PathBuf>,

    #[clap(
        long = "guide-tss-flank",
        help_heading = "Reference data",
        help = "TSS evidence search flank in bp; used only with --guide-tss",
        default_value_t = 10000,
        value_name = "BP"
    )]
    pub guide_tss_flank: u32,

    #[clap(
        long = "guide-tes-flank",
        help_heading = "Reference data",
        help = "TES evidence search flank in bp; used only with --guide-tes",
        default_value_t = 100,
        value_name = "BP"
    )]
    pub guide_tes_flank: u32,

    #[clap(
        long = "fsm-end-match-bp",
        help_heading = "Parameters",
        help = "FSM subcategory threshold for TSS/TES matching in bp; applies to FSM and ISM subcategories",
        default_value_t = 50,
        value_name = "BP"
    )]
    pub fsm_end_match_bp: i32,

    #[clap(
        long = "downstream-len",
        help = "Downstream sequence length for QC",
        help_heading = "Parameters",
        default_value_t = 20,
        value_name = "BP"
    )]
    pub downstream_len: usize,

    #[clap(
        long = "motif-search-window",
        help_heading = "Parameters",
        help = "Window size for searching polyA motifs downstream of TES in bp; applies to SQANTI3-style motif classification",
        default_value_t = 50,
        value_name = "BP"
    )]
    pub motif_search_window: usize,

    #[clap(
        long = "skip-missing-ref-chr",
        action = ArgAction::SetTrue,
        help_heading = "Auto Indexing",
        help = "Skip transcripts on seqids absent from the reference FASTA"
    )]
    pub skip_missing_ref_chr: bool,
}

// for comparison two GTFs and report F1 etc..
#[derive(Parser, Debug, Serialize, Clone)]
#[clap(
    about = "compare input GTF with base line GTF and report comparison metrics
"
)]
pub struct BenchArgs {
    #[clap(help_heading = "Input", help = "Compared GTF")]
    pub input: PathBuf,

    #[clap(
        short = 'o',
        long = "out",
        help_heading = "Output",
        help = "Output prefix"
    )]
    pub out: PathBuf,

    #[clap(
        short = 'b',
        long = "base-gtf",
        help_heading = "Baseline",
        help = "Baseline GTF"
    )]
    pub ref_gtf: PathBuf,

    #[clap(
        short = 'r',
        long = "ref-fa",
        help_heading = "Auto Indexing",
        help = "Reference FASTA"
    )]
    pub ref_fa: PathBuf,

    #[clap(
        long = "skip-missing-ref-chr",
        action = ArgAction::SetTrue,
        help_heading = "Auto Indexing",
        help = "Skip transcripts on seqids absent from the reference FASTA"
    )]
    pub skip_missing_ref_chr: bool,
}

#[derive(Debug, Clone, Serialize, ValueEnum)]
pub enum ChopMode {
    All,
    Isomatch,
}
#[derive(Parser, Debug, Serialize, Clone)]
#[clap(about = "Remove attributes from the GTF file.
")]
pub struct ChopArgs {
    #[clap(help_heading = "Input", help = "Input GTF")]
    pub input: PathBuf,

    #[clap(
        short = 'o',
        long = "out",
        help_heading = "Output",
        help = "Output prefix"
    )]
    pub out: PathBuf,

    #[clap(
        short = 'm',
        long = "mode",
        help_heading = "Parameters",
        help = "Chop mode",
        value_enum,
        default_value_t = ChopMode::Isomatch
    )]
    pub chop_mode: ChopMode,

    #[clap(
        short = 'k',
        long = "keep",
        help_heading = "Parameters",
        help = "Keep attributes"
    )]
    pub keep_attrs: Option<String>,

    #[clap(
        short = 'c',
        long = "keep-check-case",
        help_heading = "Parameters",
        help = "check letter case when match attributes",
        action = ArgAction::SetTrue,
    )]
    pub keep_check_case: bool,
}

#[derive(Subcommand, Debug, Clone, Serialize)]
pub enum ToolssArgs {
    Chop(ChopArgs),
}

fn main() {
    // set env logger level to info by default, can be overridden by RUST_LOG env var
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    match Cli::parse() {
        Cli {
            command: Commands::Index(mut args),
        } => {
            if let Err(e) = run_index(&mut args) {
                error!("{}", e);
                std::process::exit(1);
            }
        }
        Cli {
            command: Commands::Merge(args),
        } => {
            if let Err(e) = run_merge(args) {
                error!("{}", e);
                std::process::exit(1);
            }
        }
        Cli {
            command: Commands::Classify(args),
        } => {
            if let Err(e) = run_classify(args) {
                error!("{}", e);
                std::process::exit(1);
            }
        }
        Cli {
            command: Commands::Tools(args),
        } => match args {
            ToolssArgs::Chop(args) => {
                if let Err(e) = run_chop(&args) {
                    error!("{}", e);
                    std::process::exit(1);
                }
            }
        },
    }
}
