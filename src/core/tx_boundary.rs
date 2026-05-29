use std::cmp::Ordering;

use crate::core::tx_strand::ISOMSTRAND;

/// Encodes a transcript boundary as a single u64:
/// [left: 32bit | right: 30bit | strand: 2bit]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(transparent)]
pub struct TxBoundary(u64);

impl TxBoundary {
    const STRAND_MASK: u64 = 0b11;
    const RIGHT_MASK: u64 = 0x3FFF_FFFF;

    #[inline(always)]
    pub fn new(left: u32, right: u32, strand: ISOMSTRAND) -> Self {
        assert!(
            right <= Self::RIGHT_MASK as u32,
            "Right boundary must fit in 30 bits"
        );
        Self(((left as u64) << 32) | ((right as u64) << 2) | (strand.to_bit() as u64))
    }

    #[inline(always)]
    pub fn left(self) -> u32 {
        (self.0 >> 32) as u32
    }

    #[inline(always)]
    pub fn right(self) -> u32 {
        ((self.0 >> 2) & Self::RIGHT_MASK) as u32
    }

    #[inline(always)]
    pub fn strand(self) -> ISOMSTRAND {
        ISOMSTRAND::from_bit((self.0 & Self::STRAND_MASK) as u8).unwrap()
    }

    #[inline(always)]
    pub fn raw(self) -> u64 {
        self.0
    }

    /// Closed-interval overlap check: [l1, r1] overlap [l2, r2]
    ///  l1 =========== r1
    ///        l2========== r2
    ///  l2 =========== r2
    ///       l1========== r1
    /// l1 <= r2 && l2 <= r1
    #[inline(always)]
    pub fn overlaps(self, other: Self) -> bool {
        let l1 = self.0 >> 32;
        let r1 = (self.0 >> 2) & Self::RIGHT_MASK;
        let l2 = other.0 >> 32;
        let r2 = (other.0 >> 2) & Self::RIGHT_MASK;

        l1 <= r2 && l2 <= r1
    }

    /// Strand-aware closed-interval overlap
    #[inline(always)]
    pub fn overlaps_stranded(self, other: Self) -> bool {
        let l1 = self.0 >> 32;
        let r1 = (self.0 >> 2) & Self::RIGHT_MASK;
        let l2 = other.0 >> 32;
        let r2 = (other.0 >> 2) & Self::RIGHT_MASK;

        self.strand() == other.strand() && l1 <= r2 && l2 <= r1
    }

    /// 检查 self 是否完全包含 other: l1 <= l2 && r2 <= r1
    #[inline(always)]
    pub fn contains(self, other: Self) -> bool {
        let l1 = self.0 >> 32;
        let r1 = (self.0 >> 2) & Self::RIGHT_MASK;
        let l2 = other.0 >> 32;
        let r2 = (other.0 >> 2) & Self::RIGHT_MASK;

        // a <= b <==> (b.wrapping_sub(a)) >> 63 == 0
        let c1 = l2.wrapping_sub(l1) >> 63; // 0 if l1 <= l2
        let c2 = r1.wrapping_sub(r2) >> 63; // 0 if r2 <= r1
        (c1 | c2) == 0
    }
}

impl Ord for TxBoundary {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }
}

impl PartialOrd for TxBoundary {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl std::fmt::Display for TxBoundary {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let strand_ch = match self.strand() {
            ISOMSTRAND::Plus => '+',
            ISOMSTRAND::Minus => '-',
            ISOMSTRAND::Unknown => '.',
        };
        write!(f, "[{}, {}]{}", self.left(), self.right(), strand_ch)
    }
}
