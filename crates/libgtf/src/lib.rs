//! GTF parsing and isomatch index I/O.
//!
//! The crate reads `transcript` and `exon` records, preserves their 1-based
//! inclusive genomic coordinates, and emits transcripts in deterministic
//! chromosome/coordinate order. [`build_index`] writes the versioned `.isomx`
//! index and `.isoms` attribute sidecar consumed by [`index::IndexReader`].
//!
//! Library operations return structured errors and never initialize a logger
//! or terminate the process. Applications may initialize any `log` backend.

mod build;
pub mod error;
pub mod fasta;
pub mod gtf;
pub mod index;
mod io;
#[cfg(test)]
mod tests;
mod traits;
mod utils;

pub use build::{BuildConfig, BuildEvent, BuildReport, build_index, build_index_with_events};
