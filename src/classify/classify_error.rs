use crate::core::core_error::TxBaseError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum ClassifyError {
    #[error("can not k-way merge")]
    BadKWayMerge,

    #[error("TxType is not correct:{reason}")]
    TxType { reason: String },

    #[error(transparent)]
    Io(#[from] std::io::Error),

    #[error(transparent)]
    Core(#[from] TxBaseError),

    #[error("Strand is not correct:{reason}")]
    InvaidStrand { reason: String },

    #[error("Failed parse gtf: {reason}")]
    FailedParseGTF { reason: String },

}
