/// This module stores per-transcript `ISOM_SRC` sidecar data for merged GTFs.
///
/// The main `.isomx` index stores transcript structure. This sidecar stores:
/// 1. one fixed-size source-record table,
/// 2. one raw string table for source transcript IDs,
/// 3. one raw exon-diff table,
/// 4. one transcript-id -> source-record span table.
///
/// We intentionally keep the sidecar simple and append-only because the write
/// path is linear during `isomatch index`.
use std::{
    fs::File,
    io::{Seek, Write},
    path::{Path, PathBuf},
};

use crate::{
    classify::classify_error::ClassifyError,
    core::{string_pool::StringSpan, tx_type::TxType},
};

const ISOM_SRC_RECORD_DISK_SIZE: u64 = 44;

#[derive(Debug, Clone, Copy)]
pub struct IsomSrcRecord {
    pub source_file_id: u32,
    pub isom_tx_id_hash: u64, // stable hash of source transcript_id
    pub tx_type: TxType,
    pub start: u32,
    pub end: u32,
    pub donor_diff: u16,
    pub acceptor_diff: u16,
    pub src_tx_id_string: StringSpan, // source transcript_id string span
    pub exon_diffs_offset: u32,       // offset in number of ExonDiff entries
    pub exon_diffs_count: u16,
}

impl IsomSrcRecord {
    fn encode_to<W: Write>(&self, writer: &mut W) -> Result<(), ClassifyError> {
        writer.write_all(&self.source_file_id.to_le_bytes())?;
        writer.write_all(&self.isom_tx_id_hash.to_le_bytes())?;
        writer.write_all(&[tx_type_to_u8(self.tx_type)])?;
        writer.write_all(&[0u8; 3])?;
        writer.write_all(&self.start.to_le_bytes())?;
        writer.write_all(&self.end.to_le_bytes())?;
        writer.write_all(&self.donor_diff.to_le_bytes())?;
        writer.write_all(&self.acceptor_diff.to_le_bytes())?;
        writer.write_all(&self.src_tx_id_string.offset.to_le_bytes())?;
        writer.write_all(&self.src_tx_id_string.byte_len.to_le_bytes())?;
        writer.write_all(&self.exon_diffs_offset.to_le_bytes())?;
        writer.write_all(&self.exon_diffs_count.to_le_bytes())?;
        writer.write_all(&[0u8; 2])?;
        Ok(())
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ExonDiff(u64);

impl ExonDiff {
    fn pack(exon_id: u32, left_diff: i32, right_diff: i32) -> u64 {
        let mask24 = 0xFFFFFF;
        (exon_id as u64 & 0xFFFF) << 48
            | (left_diff as u64 & mask24) << 24
            | (right_diff as u64 & mask24)
    }
}

// Parse "(n exon,start diff,end diff),(n exon,start diff,end diff),..."
// into packed ExonDiff entries.
fn parse_exon_diffs(s: &str) -> Result<Option<Vec<ExonDiff>>, ClassifyError> {
    if s == "no_diff" {
        return Ok(None);
    }

    let mut result = Vec::new();
    for raw in s.split("),(") {
        let trimmed = raw.trim_matches(|c| c == '(' || c == ')');
        let nums: Vec<&str> = trimmed.splitn(3, ',').collect();
        if nums.len() != 3 {
            return Err(ClassifyError::ParseSrcRecord {
                reason: format!("invalid exon_diff tuple: {:?}", raw),
            });
        }

        let exon_id: u32 = nums[0].parse().map_err(|_| ClassifyError::ParseSrcRecord {
            reason: format!("invalid exon_id in diff: {}", nums[0]),
        })?;
        let left_diff: i32 = nums[1].parse().map_err(|_| ClassifyError::ParseSrcRecord {
            reason: format!("invalid left_diff: {}", nums[1]),
        })?;
        let right_diff: i32 = nums[2].parse().map_err(|_| ClassifyError::ParseSrcRecord {
            reason: format!("invalid right_diff: {}", nums[2]),
        })?;
        result.push(ExonDiff(ExonDiff::pack(exon_id, left_diff, right_diff)));
    }

    Ok(Some(result))
}

#[derive(Debug, Clone, Copy, Default)]
pub struct IosmSrcSpan {
    pub offset: u32, // record index in the source-record table
    pub count: u32,  // number of contiguous source records
}

pub struct IsomAttrCache {
    span_file_name: PathBuf,
    src_file: File,
    string_pool_file: File,
    exon_diffs_file: File,
    // The merged transcript ID written in the GTF (`ISOMT_n`) is parsed from
    // transcript attributes before we know whether the transcript survives
    // reference filtering. We keep the original mapping first, then project it
    // onto the dense tx_idx used by `TxBase`.
    original_tx_id_to_span: Vec<IosmSrcSpan>,
    tx_id_to_span: Vec<IosmSrcSpan>,
}

impl IsomAttrCache {
    pub fn init<P: AsRef<Path>>(path: P) -> Self {
        let base_path = path.as_ref().to_path_buf();
        let span_file_name = sidecar_path(&base_path, "isomsrc.spans");
        let src_file_name = sidecar_path(&base_path, "isomsrc.records");
        let string_pool_file_name = sidecar_path(&base_path, "isomsrc.strings");
        let exon_diffs_file_name = sidecar_path(&base_path, "isomsrc.exon_diffs");

        Self {
            span_file_name,
            src_file: File::create(&src_file_name).unwrap_or_else(|e| {
                panic!(
                    "failed to create ISOM_SRC record sidecar {}: {}",
                    src_file_name.display(),
                    e
                )
            }),
            string_pool_file: File::create(&string_pool_file_name).unwrap_or_else(|e| {
                panic!(
                    "failed to create ISOM_SRC string sidecar {}: {}",
                    string_pool_file_name.display(),
                    e
                )
            }),
            exon_diffs_file: File::create(&exon_diffs_file_name).unwrap_or_else(|e| {
                panic!(
                    "failed to create ISOM_SRC exon-diff sidecar {}: {}",
                    exon_diffs_file_name.display(),
                    e
                )
            }),
            original_tx_id_to_span: Vec::new(),
            tx_id_to_span: Vec::new(),
        }
    }

    /// Parse and append one source record from an ISOM_SRC entry
    /// (`S{file_id}:tx_id:start:end:tx_type:donor_diff:acceptor_diff:exon_diffs`).
    ///
    /// Returns the zero-based record index inside the source-record table.
    pub fn add_iso_src_from_str(&mut self, src_str: &str) -> Result<usize, ClassifyError> {
        let parts: Vec<&str> = src_str.splitn(8, ':').collect();
        if parts.len() != 8 {
            return Err(ClassifyError::ParseSrcRecord {
                reason: format!("invalid ISOM_SRC item: {}", src_str),
            });
        }

        let source_file_id = parts[0]
            .strip_prefix('S')
            .ok_or_else(|| ClassifyError::ParseSrcRecord {
                reason: format!("invalid source file id prefix: {}", parts[0]),
            })?
            .parse::<u32>()?;
        if source_file_id == 0 {
            return Err(ClassifyError::ParseSrcRecord {
                reason: "source file ids are 1-based in ISOM_SRC".to_string(),
            });
        }

        let src_tx_id = parts[1];
        let start = parts[2].parse::<u32>()?;
        let end = parts[3].parse::<u32>()?;
        let tx_type = TxType::from_str(parts[4]).ok_or_else(|| ClassifyError::TxType {
            reason: format!("unknown tx_type in ISOM_SRC: {}", parts[4]),
        })?;

        let donor_diff_raw = parts[5].parse::<u32>()?;
        let donor_diff =
            u16::try_from(donor_diff_raw).map_err(|_| ClassifyError::ParseSrcRecord {
                reason: format!("donor_diff exceeds u16: {}", donor_diff_raw),
            })?;

        let acceptor_diff_raw = parts[6].parse::<u32>()?;
        let acceptor_diff =
            u16::try_from(acceptor_diff_raw).map_err(|_| ClassifyError::ParseSrcRecord {
                reason: format!("acceptor_diff exceeds u16: {}", acceptor_diff_raw),
            })?;

        let string_offset = self.string_pool_file.stream_position()?;
        self.string_pool_file.write_all(src_tx_id.as_bytes())?;
        let src_tx_id_string = StringSpan {
            offset: u32::try_from(string_offset).map_err(|_| ClassifyError::ParseSrcRecord {
                reason: "string sidecar offset exceeded u32".to_string(),
            })?,
            byte_len: u32::try_from(src_tx_id.len()).map_err(|_| {
                ClassifyError::ParseSrcRecord {
                    reason: "source transcript id length exceeded u32".to_string(),
                }
            })?,
        };

        let exon_diffs_offset_bytes = self.exon_diffs_file.stream_position()?;
        let exon_diffs_offset = u32::try_from(
            exon_diffs_offset_bytes / std::mem::size_of::<u64>() as u64,
        )
        .map_err(|_| ClassifyError::ParseSrcRecord {
            reason: "exon-diff sidecar offset exceeded u32".to_string(),
        })?;

        let exon_diffs = parse_exon_diffs(parts[7])?;
        let exon_diffs_count = if let Some(exon_diffs) = exon_diffs {
            let exon_diffs_count =
                u16::try_from(exon_diffs.len()).map_err(|_| ClassifyError::ParseSrcRecord {
                    reason: "too many exon-diff items in one ISOM_SRC entry".to_string(),
                })?;
            for exon_diff in exon_diffs {
                self.exon_diffs_file.write_all(&exon_diff.0.to_le_bytes())?;
            }
            exon_diffs_count
        } else {
            0
        };

        let record_offset_bytes = self.src_file.stream_position()?;
        let record_index = usize::try_from(record_offset_bytes / ISOM_SRC_RECORD_DISK_SIZE)
            .map_err(|_| ClassifyError::ParseSrcRecord {
                reason: "source-record sidecar offset exceeded usize".to_string(),
            })?;

        let record = IsomSrcRecord {
            source_file_id: source_file_id - 1,
            isom_tx_id_hash: xxhash_rust::xxh3::xxh3_64(src_tx_id.as_bytes()),
            tx_type,
            start,
            end,
            donor_diff,
            acceptor_diff,
            src_tx_id_string,
            exon_diffs_offset,
            exon_diffs_count,
        };
        record.encode_to(&mut self.src_file)?;

        Ok(record_index)
    }

    /// Parse a full `ISOM_SRC` attribute string and store the contiguous span
    /// keyed by the original merged transcript id from the GTF.
    pub fn dump_isom_src_string(&mut self, s: String, tx_id: u32) -> Result<(), ClassifyError> {
        let span_start_bytes = self.src_file.stream_position()?;
        let span_start =
            u32::try_from(span_start_bytes / ISOM_SRC_RECORD_DISK_SIZE).map_err(|_| {
                ClassifyError::ParseSrcRecord {
                    reason: "source-record span offset exceeded u32".to_string(),
                }
            })?;

        let mut src_count = 0u32;
        for src_str in s.split('|').filter(|item| !item.is_empty()) {
            self.add_iso_src_from_str(src_str)?;
            src_count += 1;
        }

        ensure_span_len(&mut self.original_tx_id_to_span, tx_id as usize + 1);
        self.original_tx_id_to_span[tx_id as usize] = IosmSrcSpan {
            offset: span_start,
            count: src_count,
        };

        Ok(())
    }

    /// Project the span keyed by the merged transcript id (`ISOMT_n`) onto the
    /// dense `tx_idx` used by the indexed `TxBase`.
    pub fn project_tx_id(&mut self, original_tx_id: u32, indexed_tx_id: u32) {
        ensure_span_len(&mut self.tx_id_to_span, indexed_tx_id as usize + 1);
        self.tx_id_to_span[indexed_tx_id as usize] = self
            .original_tx_id_to_span
            .get(original_tx_id as usize)
            .copied()
            .unwrap_or_default();
    }

    pub fn finalize(mut self) -> Result<(), ClassifyError> {
        if self.tx_id_to_span.is_empty() && !self.original_tx_id_to_span.is_empty() {
            // Fallback for callers that do not need tx-id compaction.
            self.tx_id_to_span = self.original_tx_id_to_span.clone();
        }

        let mut span_file = File::create(&self.span_file_name)?;
        for span in &self.tx_id_to_span {
            span_file.write_all(&span.offset.to_le_bytes())?;
            span_file.write_all(&span.count.to_le_bytes())?;
        }

        self.src_file.flush()?;
        self.string_pool_file.flush()?;
        self.exon_diffs_file.flush()?;
        span_file.flush()?;
        Ok(())
    }
}

fn sidecar_path(base_path: &Path, suffix: &str) -> PathBuf {
    let mut path = base_path.to_path_buf();
    path.add_extension(suffix);
    path
}

fn ensure_span_len(spans: &mut Vec<IosmSrcSpan>, min_len: usize) {
    if spans.len() < min_len {
        spans.resize(min_len, IosmSrcSpan::default());
    }
}

fn tx_type_to_u8(tx_type: TxType) -> u8 {
    match tx_type {
        TxType::MONO => 0,
        TxType::ALLC => 1,
        TxType::PRTC => 2,
        TxType::NOTC => 3,
    }
}
