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
}
