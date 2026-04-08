use crate::core::tx_base_error::TxBaseError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum MergeError {
    #[error("can not k-way merge")]
    BadKWayMerge,

    #[error("TxType is not correct:{reason}")]
    TxType { reason: String },

    #[error(transparent)]
    Io(#[from] std::io::Error),

    #[error(transparent)]
    Core(#[from] TxBaseError),
}
