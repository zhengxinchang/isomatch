use std::{
    fs::File,
    io::{BufWriter, Seek, SeekFrom, Write},

};


use crate::{ index::index_error::IndexError};

use libgtf::index::ISOMS_VERSION;
/// Sidecar file layout (.isomattr):
///   [Header:     41 bytes - magic(5) + version(4) + md5(16) + total_tx_n(8) + span_table_off(8)]
///   [Blob:       variable  — per-tx zstd-compressed attr bytes, written in tx_gidx order]
///   [Span table: N x 12B  - RawStringSpan entries indexed by tx_gidx]

#[derive(Clone, Copy, Default)]
pub struct RawStringSpan {
    pub offset: u64, // absolute byte offset of the compressed blob in the file
    pub length: u32, // byte length of the compressed blob
}

const MAGIC: [u8; 5] = *b"ISOMS";

// magic(5) + version(4) + md5(16) + total_tx_n(8) + span_table_off(8) = 41
const HEADER_SIZE: usize = 41;

pub struct AttrIndexBuilder {
    magic: [u8; 5],
    version: u32,
    md5: [u8; 16],
    total_tx_n: u64,
    current_tx_n: usize,
    blob_offset: usize, // current write cursor; starts at HEADER_SIZE, grows with each dump_attr
    tx_idx_to_raw_attr_span: Vec<RawStringSpan>,
    file: BufWriter<File>,
}

impl AttrIndexBuilder {
    pub(crate) fn new(file: File, total_tx_n: u64, md5: &[u8; 16]) -> Result<Self, IndexError> {
        let mut file = BufWriter::new(file);
        // Reserve space for header; will be overwritten in finish()
        file.write_all(&[0u8; HEADER_SIZE])
            .map_err(|e| IndexError::FailReadIndex {
                reason: e.to_string(),
            })?;
        Ok(Self {
            magic: MAGIC,
            version: ISOMS_VERSION,
            md5: *md5,
            total_tx_n,
            current_tx_n: 0,
            blob_offset: HEADER_SIZE,
            tx_idx_to_raw_attr_span: vec![
                RawStringSpan::default();
                usize::try_from(total_tx_n).map_err(|_| {
                    IndexError::FailReadIndex {
                        reason: format!("total_tx_n {} exceeded usize", total_tx_n),
                    }
                })?
            ],
            file,
        })
    }

    /// Compress `data` with zstd and append to the blob section.
    /// Records the (offset, length) span for `tx_gidx` in the span table.
    pub fn dump_attr(&mut self, data: Vec<u8>, tx_gidx: u64) -> Result<usize, IndexError> {
        let idx = usize::try_from(tx_gidx).map_err(|_| IndexError::FailReadIndex {
            reason: format!("tx_gidx {} exceeded usize", tx_gidx),
        })?;
        if tx_gidx >= self.total_tx_n {
            return Err(IndexError::FailReadIndex {
                reason: format!(
                    "tx_gidx {} out of range (total {})",
                    tx_gidx, self.total_tx_n
                ),
            });
        }
        let compressed =
            zstd::encode_all(data.as_slice(), 3).map_err(|e| IndexError::FailReadIndex {
                reason: e.to_string(),
            })?;
        let offset = u64::try_from(self.blob_offset).map_err(|_| IndexError::FailReadIndex {
            reason: format!("blob offset {} exceeded u64", self.blob_offset),
        })?;
        let length = u32::try_from(compressed.len()).map_err(|_| IndexError::FailReadIndex {
            reason: format!("compressed attr length {} exceeded u32", compressed.len()),
        })?;
        self.file
            .write_all(&compressed)
            .map_err(|e| IndexError::FailReadIndex {
                reason: e.to_string(),
            })?;
        self.tx_idx_to_raw_attr_span[idx] = RawStringSpan { offset, length };
        self.blob_offset += compressed.len();
        self.current_tx_n += 1;
        Ok(compressed.len())
    }

    /// Finalize: seek back to write the real header, then append the span table.
    pub fn finish(mut self) -> Result<(), IndexError> {
        // Use a fn-pointer so it is Copy and can be reused across map_err calls.
        let e: fn(std::io::Error) -> IndexError = |err| IndexError::FailReadIndex {
            reason: err.to_string(),
        };

        let span_table_off = self.blob_offset as u64;

        // Overwrite the placeholder header at position 0.
        self.file.seek(SeekFrom::Start(0)).map_err(e)?;
        self.file.write_all(&self.magic).map_err(e)?;
        self.file
            .write_all(&self.version.to_le_bytes())
            .map_err(e)?;
        self.file.write_all(&self.md5).map_err(e)?;
        self.file
            .write_all(&self.total_tx_n.to_le_bytes())
            .map_err(e)?;
        self.file
            .write_all(&span_table_off.to_le_bytes())
            .map_err(e)?;

        // Append span table after the blob section.
        self.file.seek(SeekFrom::Start(span_table_off)).map_err(e)?;
        for span in &self.tx_idx_to_raw_attr_span {
            self.file.write_all(&span.offset.to_le_bytes()).map_err(e)?;
            self.file.write_all(&span.length.to_le_bytes()).map_err(e)?;
        }

        self.file.flush().map_err(e)
    }
}

