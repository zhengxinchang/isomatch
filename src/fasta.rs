use std::{
    io::{BufRead, Seek, SeekFrom},
    path::Path,
};

use noodles_core::{Position, region::Interval};
use noodles_fasta as fasta;
use rustc_hash::FxHashMap;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum FastaError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("duplicate seqid in FASTA index: '{0}'")]
    DuplicateSeqId(String),
    #[error("seqid not found: '{0}'")]
    SeqIdNotFound(String),
    #[error("region [{start}, {end}) out of bounds for '{seqid}' (len={seq_len})")]
    OutOfBounds {
        seqid: String,
        start: usize,
        end: usize,
        seq_len: usize,
    },
    #[error("invalid position: {0}")]
    InvalidPosition(String),
}

pub enum FaType {
    Ref,
    Seq,
}

pub struct FastaReader {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<std::fs::File>>,

    /// Maps raw FASTA seqid → index into the fai::Index Vec.
    name_to_idx: FxHashMap<String, u32>,
}

impl FastaReader {
    /// Open a FASTA file. The `.fai` index must exist at `<path>.fai`.
    pub fn open<P: AsRef<Path>>(path: P, _fa_type: FaType) -> Result<Self, FastaError> {
        let reader = fasta::io::indexed_reader::Builder::default().build_from_path(&path)?;

        let mut name_to_idx: FxHashMap<String, u32> = FxHashMap::default();

        for (idx, record) in reader.index().as_ref().iter().enumerate() {
            let name = String::from_utf8_lossy(record.name()).into_owned();

            if name_to_idx.contains_key(&name) {
                return Err(FastaError::DuplicateSeqId(name));
            }

            name_to_idx.insert(
                name,
                u32::try_from(idx).expect("fai index exceeds u32::MAX entries"),
            );
        }

        Ok(Self {
            reader,
            name_to_idx,
        })
    }

    /// Fetch sequence bytes for `seqid` in `[start, end)` (0-based, half-open).
    ///
    /// Returns forward strand by default. When `neg_strand = true`, returns
    /// the reverse complement of the fetched region.
    /// strand 0 == plus; 1 == minus
    pub fn fetch(
        &mut self,
        seqid: &str,
        start: usize,
        end: usize,
        // strand: u8,
        _trim_chr: bool,
    ) -> Result<Vec<u8>, FastaError> {
        let idx =
            self.name_to_idx
                .get(seqid)
                .copied()
                .ok_or_else(|| FastaError::SeqIdNotFound(seqid.to_string()))? as usize;

        let fai_record = &self.reader.index().as_ref()[idx];
        let seq_len = fai_record.length() as usize;

        if end > seq_len {
            return Err(FastaError::OutOfBounds {
                seqid: seqid.to_string(),
                start,
                end,
                seq_len,
            });
        }

        if start == end {
            return Ok(Vec::new());
        }

        // noodles intervals are 1-based, closed [start+1, end]
        let interval_start = Position::try_from(start + 1)
            .map_err(|e| FastaError::InvalidPosition(e.to_string()))?;
        let interval_end =
            Position::try_from(end).map_err(|e| FastaError::InvalidPosition(e.to_string()))?;
        let interval = Interval::from(interval_start..=interval_end);

        // Compute the file offset directly from the fai record — O(1), no linear scan.
        let file_pos = fai_record.query(interval).map_err(FastaError::Io)?;

        let len = end - start;

        // Seek and read using the inner BufReader directly.
        self.reader
            .get_mut()
            .seek(SeekFrom::Start(file_pos))
            .map_err(FastaError::Io)?;

        let seq = read_sequence_stripped(self.reader.get_mut(), len).map_err(FastaError::Io)?;

        Ok(seq)
        // if strand == 1 {
        //     Ok(rev_comp(&seq))
        // } else {
        //     Ok(seq)
        // }
    }

    pub fn fetch_all(&mut self, seqid: &str, trim_chr: bool) -> Result<Vec<u8>, FastaError> {
        let len = self
            .seq_len(seqid)
            .ok_or_else(|| FastaError::SeqIdNotFound(seqid.to_string()))?;

        self.fetch(seqid, 0, len, trim_chr)
    }

    pub fn contains(&self, seqid: &str) -> bool {
        self.name_to_idx.contains_key(seqid)
    }

    pub fn seq_len(&self, seqid: &str) -> Option<usize> {
        let idx = *self.name_to_idx.get(seqid)? as usize;
        Some(self.reader.index().as_ref()[idx].length() as usize)
    }

    pub fn seqids(&self) -> impl Iterator<Item = &str> {
        self.name_to_idx.keys().map(String::as_str)
    }
}

/// Read exactly `len` bases from the reader, skipping newline characters.
fn read_sequence_stripped<R: BufRead>(reader: &mut R, len: usize) -> std::io::Result<Vec<u8>> {
    let mut buf = Vec::with_capacity(len);
    while buf.len() < len {
        let need = len - buf.len();
        let available = reader.fill_buf()?;
        if available.is_empty() {
            break;
        }
        let chunk = &available[..available.len().min(
            need + available
                .iter()
                .filter(|&&b| b == b'\n' || b == b'\r')
                .count(),
        )];
        for &b in chunk {
            if b != b'\n' && b != b'\r' {
                buf.push(b);
                if buf.len() == len {
                    break;
                }
            }
        }
        let consumed = chunk.len();
        reader.consume(consumed);
    }
    Ok(buf)
}
