use thiserror::Error;

#[derive(Error, Debug)]
pub enum ToolError {
    #[error(transparent)]
    Io(#[from] std::io::Error),

    #[error("Strand is not correct:{reason}")]
    InvaidStrand { reason: String },

    #[error("Failed parse gtf: {reason}")]
    FailedParseGTF { reason: String },

    #[error(transparent)]
    ParseIntError(#[from] std::num::ParseIntError),

    #[error("Can not parse isomatch merged GTF:{reason}")]
    ReadMergedGTFFailed { reason: String },

    #[error("Invalid path, can not extract file name: {path}")]
    InvalidPath { path: String },

    #[error(transparent)]
    Json(#[from] serde_json::Error),
}
