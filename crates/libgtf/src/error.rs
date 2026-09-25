use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error(transparent)]
    Io(#[from] std::io::Error),

    #[error(transparent)]
    IndexData(#[from] IndexDataError),

    #[error(transparent)]
    Fasta(#[from] FastaError),

    #[error("failed to build index: {reason}")]
    IndexBuild { reason: String },

    #[error("failed to build transcript {tx_id}: {reason}")]
    BuildTranscript { tx_id: String, reason: String },

    #[error("invalid GTF at line {line}: {reason}")]
    InvalidGtfLine { line: usize, reason: String },

    #[error("invalid GTF record: {reason}")]
    InvalidGtfRecord {
        field: Option<usize>,
        reason: String,
    },

    #[error("invalid GTF field {field}: {reason}")]
    InvalidGtfField { field: usize, reason: String },

    #[error("invalid GTF chromosome: {reason}")]
    InvalidGtfChr { reason: String },

    #[error("invalid index: {reason}")]
    InvalidIndex { reason: String },

    #[error("unsupported index version {found}; expected {expected}")]
    UnsupportedVersion { found: u32, expected: u32 },

    #[error("reference FASTA is missing sequence {seqid:?}")]
    MissingReferenceSequence { seqid: String },

    #[error("reference FASTA is missing seqids required by the GTF: {seqids:?}")]
    MissingReferenceSequences { seqids: Vec<String> },

    #[error("no indexable seqids remain after filtering against the reference FASTA")]
    NoIndexableSequences,

    #[error("invalid coordinate {start}-{end}")]
    InvalidCoordinate { start: u32, end: u32 },

    #[error("invalid stand (u8){strand}")]
    InvalidStrand { strand: u8 },

    #[error("cannot create temporary files under {path}")]
    TempDirectory {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("GTF must contain at least one transcript record")]
    MissingTranscriptRecord,

    #[error(
        "GTF transcript/exon chromosome mismatch. Transcript-only seqids: {transcript_only:?}; exon-only seqids: {exon_only:?}"
    )]
    TranscriptExonChromMismatch {
        transcript_only: Vec<String>,
        exon_only: Vec<String>,
    },

    #[error("inconsistent transcript {tx_id}: {reason}")]
    InconsistentTranscript { tx_id: String, reason: String },

    #[error("Chromosome not found: {name}")]
    ChromosomeNotFound { name: String },
}

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum IndexDataError {
    #[error("strand must be 0 (+), 1 (-), or 2 (unknown), got {strand}")]
    InvalidStrand { strand: u8 },

    #[error("tx start {start} is greater than end {end}")]
    InvalidBounds { start: u32, end: u32 },

    #[error("n_exons must be at least 1, got {n_exons}")]
    InvalidExonCount { n_exons: u16 },

    #[error("junction coordinates must be strictly increasing")]
    JunctionsNotStrictlyIncreasing,

    #[error("too many junction coordinates for one transcript: {count}")]
    TooManyJunctions { count: usize },

    #[error("junction pool is too large to address with u64 offsets")]
    PoolTooLarge,

    #[error("junction pool mismatch: pool chrom = ({pool_chrom_id}), tx chrom = ({tx_chrom_id})")]
    PoolMismatch {
        pool_chrom_id: u16,
        tx_chrom_id: u16,
        tx_strand: u8,
    },

    #[error("invalid junction span offset={offset} count={count} for pool length {pool_len}")]
    InvalidSpan {
        offset: u32,
        count: u16,
        pool_len: usize,
    },

    #[error("invalid intern id {id} not found in string pool")]
    InvalidInternId { id: u64 },

    #[error("invalid splice site: {site}")]
    InvalidSpliceSite { site: String },

    #[error("string pool exceeded u32-addressable size")]
    StringPoolTooLarge,

    #[error("IO error: {0}")]
    Io(String),

    #[error("invalid encoding: {msg}")]
    InvalidEncoding { msg: String },
}

impl IndexDataError {
    pub fn io(err: std::io::Error) -> Self {
        Self::Io(err.to_string())
    }
}

#[derive(Debug, Error)]
pub enum FastaError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("duplicate seqid in FASTA index: '{0}'")]
    DuplicateSeqId(String),
    #[error("seqid not found: '{0}'")]
    SeqIdNotFound(String),
    #[error("region [{start}, {end}) out of bounds for '{seqid}' (len={seq_len})")]
    OutOfBounds {
        seqid: String,
        start: usize,
        end: usize,
        seq_len: usize,
    },
    #[error("invalid position: {0}")]
    InvalidPosition(String),

    #[error("Failed fetch sequence: {reason}")]
    FetchSeqFailed { reason: String },
}
