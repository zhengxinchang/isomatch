use std::{
    fs::File,
    io::{BufWriter, Read, Seek, SeekFrom, Write},
    path::Path,
};

use log::error;

use crate::{constants::ISOMS_VERSION, index::index_error::IndexError};

/// Sidecar file layout (.isomattr):
///   [Header:     20 bytes  — magic(7) + version(1) + total_tx_n(4) + span_table_off(8)]
///   [Blob:       variable  — per-tx zstd-compressed attr bytes, written in tx_gidx order]
///   [Span table: N × 8B   — RawStringSpan entries indexed by tx_gidx]

#[derive(Clone, Copy, Default)]
pub struct RawStringSpan {
    pub offset: u32, // absolute byte offset of the compressed blob in the file
    pub length: u32, // byte length of the compressed blob
}

const MAGIC: [u8; 5] = *b"ISOMS";

// magic(5) + version(1) + total_tx_n(4) + span_table_off(8) = 18
const HEADER_SIZE: usize = 18;

pub struct AttrIndexBuilder {
    magic: [u8; 5],
    version: u8,
    total_tx_n: usize,
    current_tx_n: usize,
    blob_offset: usize, // current write cursor; starts at HEADER_SIZE, grows with each dump_attr
    tx_idx_to_raw_attr_span: Vec<RawStringSpan>,
    file: BufWriter<File>,
}

impl AttrIndexBuilder {
    pub fn init<P: AsRef<Path>>(path: P, total_tx_n: usize) -> Result<Self, IndexError> {
        let file = File::create(path).map_err(|e| IndexError::FailReadIndex {
            reason: e.to_string(),
        })?;
        let mut file = BufWriter::new(file);
        // Reserve space for header; will be overwritten in finish()
        file.write_all(&[0u8; HEADER_SIZE])
            .map_err(|e| IndexError::FailReadIndex {
                reason: e.to_string(),
            })?;
        Ok(Self {
            magic: MAGIC,
            version: ISOMS_VERSION,
            total_tx_n,
            current_tx_n: 0,
            blob_offset: HEADER_SIZE,
            tx_idx_to_raw_attr_span: vec![RawStringSpan::default(); total_tx_n],
            file,
        })
    }

    /// Compress `data` with zstd and append to the blob section.
    /// Records the (offset, length) span for `tx_gidx` in the span table.
    pub fn dump_attr(&mut self, data: Vec<u8>, tx_gidx: u32) -> Result<usize, IndexError> {
        let idx = tx_gidx as usize;
        if idx >= self.total_tx_n {
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
        let offset = u32::try_from(self.blob_offset).map_err(|_| IndexError::FailReadIndex {
            reason: format!("blob offset {} exceeded u32", self.blob_offset),
        })?;
        let length = compressed.len() as u32;
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
        self.file.write_all(&[self.version]).map_err(e)?;
        self.file
            .write_all(&(self.total_tx_n as u32).to_le_bytes())
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

pub struct AttrIndexHeader {
    _magic: [u8; 5],
    _version: u8,
    total_tx_n: u32,
    span_table_off: u64,
}

pub struct AttrIndexReader {
    file: File,
    header: AttrIndexHeader,
    // total_tx_n: u32,
    // span_table_off: u64,
}

impl AttrIndexReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, IndexError> {
        // let e: fn(std::io::Error) -> IndexError = |err| IndexError::FailReadIndex {
        //     reason: err.to_string(),
        // };

        let mut file = File::open(path).map_err(|e| IndexError::FailReadIndex {
            reason: "Can not read AttrIndex file".to_string(),
        })?;

        // let mut header =

        let mut magic = [0u8; 5];
        file.read_exact(&mut magic)
            .map_err(|e| IndexError::FailReadIndex {
                reason: "Can not read magic in AttrIndex file".to_string(),
            })?;
        if magic != MAGIC {
            return Err(IndexError::FailReadIndex {
                reason: format!("invalid magic: expected ISOMS, got {:?}", magic),
            });
        }

        let mut version = [0u8; 1];
        file.read_exact(&mut version)
            .map_err(|e| IndexError::FailReadIndex {
                reason: "Can not read version in AttrIndex file".to_string(),
            })?;
        if version[0] != ISOMS_VERSION {
            error!(
                "The isomx version ({}) is outdated, please rebuild the index.",
                version[0]
            );
            std::process::exit(1);
        }

        let mut buf4 = [0u8; 4];
        file.read_exact(&mut buf4)
            .map_err(|e| IndexError::FailReadIndex {
                reason: "Can not read total tx number in AttrIndex file".to_string(),
            })?;
        let total_tx_n = u32::from_le_bytes(buf4);

        let mut buf8 = [0u8; 8];
        file.read_exact(&mut buf8)
            .map_err(|e| IndexError::FailReadIndex {
                reason: "Can not read span table offset in AttrIndex file".to_string(),
            })?;
        let span_table_off = u64::from_le_bytes(buf8);

        let header = AttrIndexHeader {
            _magic: magic,
            _version: version[0],
            total_tx_n,
            span_table_off,
        };

        Ok(Self { file, header })
    }

    /// Returns the decompressed attr bytes for `tx_idx`, or `None` if not set.
    pub fn get_attr(&mut self, tx_idx: u32) -> Result<Option<Vec<u8>>, IndexError> {
        if tx_idx >= self.header.total_tx_n {
            return Ok(None);
        }
        let e: fn(std::io::Error) -> IndexError = |err| IndexError::FailReadIndex {
            reason: err.to_string(),
        };

        // Read the RawStringSpan for this tx_idx from the span table.
        let span_entry_off = self.header.span_table_off + tx_idx as u64 * 8;
        self.file.seek(SeekFrom::Start(span_entry_off)).map_err(e)?;
        let mut buf4 = [0u8; 4];
        self.file.read_exact(&mut buf4).map_err(e)?;
        let offset = u32::from_le_bytes(buf4);
        self.file.read_exact(&mut buf4).map_err(e)?;
        let length = u32::from_le_bytes(buf4);

        if length == 0 {
            return Ok(None);
        }

        // Read and decompress the blob.
        self.file.seek(SeekFrom::Start(offset as u64)).map_err(e)?;
        let mut compressed = vec![0u8; length as usize];
        self.file.read_exact(&mut compressed).map_err(e)?;
        let data =
            zstd::decode_all(compressed.as_slice()).map_err(|err| IndexError::FailReadIndex {
                reason: err.to_string(),
            })?;

        Ok(Some(data))
    }
}
