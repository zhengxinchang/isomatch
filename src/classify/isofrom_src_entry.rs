use crate::{
    classify::classify_error::ClassifyError,
    core::{
        string_pool::{StringPool, StringSpan},
        tx_type::TxType,
    },
};

struct IsomSrcRecord {
    source_file_id: u32,
    tx_type: TxType,
    start: u32,
    end: u32,
    donor_diff: u16,
    acceptor_diff: u16,
    src_tx_id: StringSpan,
    exon_diffs_offset: u32,
    exon_diffs_count: u16,
}

pub struct ExonDiff(u64);

impl ExonDiff {
    fn pack(exon_id: u32, left_diff: i32, right_diff: i32) -> u64 {
        let mask24 = 0xFFFFFF;
        (exon_id as u64 & 0xFFFF) << 48
            | (left_diff as u64 & mask24) << 24
            | (right_diff as u64 & mask24)
    }

    fn unpack(packed: u64) -> (u32, i32, i32) {
        let exon_id = (packed >> 48) as u32;

        let left_diff: i32 = ((packed >> 24) & 0xFFFFFF) as i32;
        let left_diff = (left_diff << 8) >> 8; // sign-extend 24-bit

        let right_diff: i32 = (packed & 0xFFFFFF) as i32;
        let right_diff: i32 = (right_diff << 8) >> 8;

        (exon_id, left_diff, right_diff)
    }
}

// Parse "(n,l,r),(n,l,r),..." into packed ExonDiff entries.
fn parse_exon_diffs(s: &str) -> Result<Vec<ExonDiff>, ClassifyError> {
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
    Ok(result)
}

/// Pool for all IsomSrcRecords parsed from ISOM_SRC attributes.
pub struct IsomSrcPool {
    records: Vec<IsomSrcRecord>,
    records_size: usize,
    exon_diff: Vec<ExonDiff>,
    string_pool: StringPool,
}

impl IsomSrcPool {
    pub fn init() -> Self {
        IsomSrcPool {
            records: Vec::new(),
            records_size: 0,
            exon_diff: Vec::new(),
            string_pool: StringPool::new(),
        }
    }

    /// Parse and add one source record from an ISOM_SRC entry (a single pipe-separated segment).
    ///
    /// Format: `S{file_id}:{tx_id}:{start}:{end}:{tx_type}:{donor_diff}:{acceptor_diff}:{exon_diffs}`
    /// `file_id` is 1-indexed in the GTF; stored here as 0-indexed.
    pub fn add_iso_src_from_str(&mut self, src_str: &str) -> Result<usize, ClassifyError> {
        // Split file_id from the front, then the 6 fixed trailing fields from the right.
        // This correctly handles tx_ids that contain ':'.
        let (file_id_raw, rest) =
            src_str
                .split_once(':')
                .ok_or_else(|| ClassifyError::ParseSrcRecord {
                    reason: format!("missing ':' in src record: {:?}", src_str),
                })?;

        // rsplitn(7) from the right yields 7 segments; after reversing:
        //   [0]=tx_id  [1]=start  [2]=end  [3]=tx_type  [4]=donor_diff  [5]=acceptor_diff  [6]=exon_diffs
        let mut rparts: Vec<&str> = rest.rsplitn(7, ':').collect();
        if rparts.len() != 7 {
            return Err(ClassifyError::ParseSrcRecord {
                reason: format!(
                    "expected 7 right-split fields, got {}: {:?}",
                    rparts.len(),
                    src_str
                ),
            });
        }
        rparts.reverse();

        // "S1" → 0-indexed file id
        if !file_id_raw.starts_with('S') {
            return Err(ClassifyError::ParseSrcRecord {
                reason: format!("file_id field must start with 'S', got: {}", file_id_raw),
            });
        }
        let raw_id: u32 = file_id_raw[1..]
            .parse()
            .map_err(|_| ClassifyError::ParseSrcRecord {
                reason: format!("invalid file_id: {}", file_id_raw),
            })?;
        let source_file_id = raw_id.saturating_sub(1);

        let src_tx_id_span = self
            .string_pool
            .add(rparts[0])
            .map_err(ClassifyError::Core)?;

        let start: u32 = rparts[1]
            .parse()
            .map_err(|_| ClassifyError::ParseSrcRecord {
                reason: format!("invalid start: {}", rparts[1]),
            })?;

        let end: u32 = rparts[2]
            .parse()
            .map_err(|_| ClassifyError::ParseSrcRecord {
                reason: format!("invalid end: {}", rparts[2]),
            })?;

        let tx_type = TxType::from_str(rparts[3]).ok_or_else(|| ClassifyError::TxType {
            reason: format!("unknown tx_type: {}", rparts[3]),
        })?;

        let donor_diff: u16 = rparts[4]
            .parse::<u32>()
            .map_err(|_| ClassifyError::ParseSrcRecord {
                reason: format!("invalid donor_diff: {}", rparts[4]),
            })
            .map(|v| v.min(u16::MAX as u32) as u16)?;

        let acceptor_diff: u16 = rparts[5]
            .parse::<u32>()
            .map_err(|_| ClassifyError::ParseSrcRecord {
                reason: format!("invalid acceptor_diff: {}", rparts[5]),
            })
            .map(|v| v.min(u16::MAX as u32) as u16)?;

        let exon_diffs = if rparts[6] == "no_diff" {
            Vec::new()
        } else {
            parse_exon_diffs(rparts[6])?
        };

        let record = IsomSrcRecord {
            source_file_id,
            tx_type,
            start,
            end,
            donor_diff,
            acceptor_diff,
            src_tx_id: src_tx_id_span,
            exon_diffs_offset: self.exon_diff.len() as u32,
            exon_diffs_count: exon_diffs.len() as u16,
        };

        self.exon_diff.extend(exon_diffs);
        self.records.push(record);
        self.records_size += 1;
        Ok(self.records_size - 1)
    }

    pub fn add_iso_src_str_vec(
        &mut self,
        src_strs: Vec<&str>,
    ) -> Result<IsoSrcSpan, ClassifyError> {
        let current_idx = self.records_size;
        for src_str in &src_strs {
            let _ = self.add_iso_src_from_str(src_str)?;
        }
        Ok(IsoSrcSpan {
            offset: current_idx,
            count: src_strs.len(),
        })
    }
}

pub struct IsoSrcSpan {
    offset: usize,
    count: usize,
}

// -- QPTIR
//     TxBase
//     IsoSrcSpan

// -- RPTIR
//     TxBase
//     OtherInfo
//         Reference IDs

// IsoSrcPool
