use std::convert::Infallible;

use thiserror::Error;

use crate::utils::UtilsError;

#[derive(Error, Debug)]
pub enum ToolError {
    #[error(transparent)]
    Io(#[from] std::io::Error),

    #[error(transparent)]
    Utils(#[from] UtilsError),

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

    #[error(transparent)]
    ParseStr(#[from] Infallible),
}
