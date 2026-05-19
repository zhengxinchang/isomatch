use crate::core::{string_pool::StringSpan, tx_type::TxType};

struct IsomSrcRecord {
    source_file_id: u8,
    tx_type: TxType, // TxType enum packed
    start: u32,
    end: u32,
    donor_diff: u32,
    acceptor_diff: u32,
    tx_id_span: StringSpan, // 复用已有 StringPool
    exon_diffs_offset: u32,
    exon_diffs_count: u16,
}

impl IsomSrcRecord {

    pub fn from_str(s:&str) -> Self {

        
    }
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
        let left_diff = (left_diff << 8) >> 8; // 24-bit 符号扩展

        let right_diff: i32 = (packed & 0xFFFFFF) as i32;
        let right_diff: i32 = (right_diff << 8) >> 8;

        (exon_id, left_diff, right_diff)
    }
}


pub struct IsomSrcPool {
    records:Vec<IsomSrcRecord>,
    exon_diff:Vec<ExonDiff>
}