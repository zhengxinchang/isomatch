use thiserror::Error;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum TxBaseError {
    #[error("strand must be 0 (+) or 1 (-), got {strand}")]
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
    InvalidInternId { id: u32 },

    #[error("string pool exceeded u32::MAX entries")]
    StringPoolTooLarge,

    #[error("IO error: {0}")]
    Io(String),

    #[error("invalid encoding: {msg}")]
    InvalidEncoding { msg: String },
}

impl TxBaseError {
    pub fn io(err: std::io::Error) -> Self {
        Self::Io(err.to_string())
    }
}
