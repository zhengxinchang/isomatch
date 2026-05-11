use thiserror::Error;

use crate::fasta::FastaError;

#[derive(Error, Debug)]
pub enum IndexError {
    #[error("Failed to add GTFTx {id} to ChromBlockBuilder. Reason: {reason}")]
    AddGTFTx { id: String, reason: String },

    #[error("Failed to add GTFTx {id} to JunctionPool. Reason: {reason}")]
    JunctionPoolAdd { id: String, reason: String },

    #[error("Failed to add GTFTx {id} to StringPool. Reason: {reason}")]
    StringPoolAdd { id: String, reason: String },

    #[error("Failed to fetch sequence from the Fasta file. Reason: {reason}")]
    FetchSeqFailed { reason: String },

    #[error(transparent)]
    Fasta(#[from] FastaError),
}
