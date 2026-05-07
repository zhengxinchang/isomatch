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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roundtrip() {
        let b = TxBoundary::new(1000, 2000, ISOMSTRAND::Plus);
        assert_eq!(b.left(), 1000);
        assert_eq!(b.right(), 2000);
        assert_eq!(b.strand(), ISOMSTRAND::Plus);
    }

    #[test]
    fn test_overlap() {
        let a = TxBoundary::new(100, 200, ISOMSTRAND::Plus);
        let b = TxBoundary::new(150, 300, ISOMSTRAND::Plus);
        let c = TxBoundary::new(200, 400, ISOMSTRAND::Plus); // [200, 400] 与 [100, 200] overlap（闭区间）
        let d = TxBoundary::new(50, 80, ISOMSTRAND::Plus);

        assert!(a.overlaps(b)); // [100,200] ∩ [150,300] ✓
        assert!(b.overlaps(a)); // 对称
        assert!(a.overlaps(c)); // [100,200] ∩ [200,400] ✓
        assert!(!a.overlaps(d)); // [100,200] ∩ [50,80]   ✗
    }

    #[test]
    fn test_contains() {
        let outer = TxBoundary::new(100, 500, ISOMSTRAND::Plus);
        let inner = TxBoundary::new(200, 400, ISOMSTRAND::Plus);
        let partial = TxBoundary::new(200, 600, ISOMSTRAND::Plus);

        assert!(outer.contains(inner));
        assert!(!outer.contains(partial));
        assert!(!inner.contains(outer));
    }

    #[test]
    fn test_sort_by_left() {
        let mut v = vec![
            TxBoundary::new(500, 600, ISOMSTRAND::Plus),
            TxBoundary::new(100, 200, ISOMSTRAND::Plus),
            TxBoundary::new(100, 150, ISOMSTRAND::Plus),
        ];
        v.sort();
        assert_eq!(v[0].left(), 100);
        assert_eq!(v[0].right(), 150); // 同 left 时 right 小的在前
        assert_eq!(v[1].right(), 200);
        assert_eq!(v[2].left(), 500);
    }
}
