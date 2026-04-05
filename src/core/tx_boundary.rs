use std::cmp::Ordering;

/// Encodes a transcript boundary as a single u64:
/// [left: 32bit | right: 31bit | strand: 1bit]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(transparent)]
pub struct TxBoundary(u64);

impl TxBoundary {
    /// Strand encoding
    pub const PLUS: u8 = 0;
    pub const MINUS: u8 = 1;

    #[inline(always)]
    pub fn new(left: u32, right: u32, strand: u8) -> Self {
        assert!(strand <= 1, "Strand must be 0 (plus) or 1 (minus)");
        assert!(right <= 0x7FFF_FFFF, "Right boundary must fit in 31 bits");
        Self(((left as u64) << 32) | ((right as u64) << 1) | (strand as u64))
    }

    #[inline(always)]
    pub fn left(self) -> u32 {
        (self.0 >> 32) as u32
    }

    #[inline(always)]
    pub fn right(self) -> u32 {
        ((self.0 >> 1) & 0x7FFF_FFFF) as u32
    }

    #[inline(always)]
    pub fn strand(self) -> u8 {
        (self.0 & 1) as u8
    }

    #[inline(always)]
    pub fn raw(self) -> u64 {
        self.0
    }

    /// Branchless overlap check: [l1, r1) overlap [l2, r2)  
    ///  l1 =========== r1
    ///        l2========== r2
    ///  l2 =========== r2
    ///       l1========== r1
    /// l1 < r2 && l2 < r1
    #[inline(always)]
    pub fn overlaps(self, other: Self) -> bool {
        let l1 = self.0 >> 32;
        let r1 = (self.0 >> 1) & 0x7FFF_FFFF;
        let l2 = other.0 >> 32;
        let r2 = (other.0 >> 1) & 0x7FFF_FFFF;

        // wrapping_sub + 符号位: branchless 比较
        // 对于非负整数 a < b <==> (a.wrapping_sub(b)) >> 63 == 1
        let c1 = l1.wrapping_sub(r2) >> 63;
        let c2 = l2.wrapping_sub(r1) >> 63;
        (c1 & c2) != 0
    }

    /// Strand-aware overlap
    #[inline(always)]
    pub fn overlaps_stranded(self, other: Self) -> bool {
        let l1 = self.0 >> 32;
        let r1 = (self.0 >> 1) & 0x7FFF_FFFF;
        let l2 = other.0 >> 32;
        let r2 = (other.0 >> 1) & 0x7FFF_FFFF;

        let c1 = l1.wrapping_sub(r2) >> 63;
        let c2 = l2.wrapping_sub(r1) >> 63;
        let strand_match = 1 - ((self.0 ^ other.0) & 1);

        (c1 & c2 & strand_match) != 0
    }

    /// 检查 self 是否完全包含 other: l1 <= l2 && r2 <= r1
    #[inline(always)]
    pub fn contains(self, other: Self) -> bool {
        let l1 = self.0 >> 32;
        let r1 = (self.0 >> 1) & 0x7FFF_FFFF;
        let l2 = other.0 >> 32;
        let r2 = (other.0 >> 1) & 0x7FFF_FFFF;

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
        let strand_ch = if self.strand() == Self::PLUS {
            '+'
        } else {
            '-'
        };
        write!(f, "[{}, {}){}", self.left(), self.right(), strand_ch)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roundtrip() {
        let b = TxBoundary::new(1000, 2000, TxBoundary::PLUS);
        assert_eq!(b.left(), 1000);
        assert_eq!(b.right(), 2000);
        assert_eq!(b.strand(), 0);
    }

    #[test]
    fn test_overlap() {
        let a = TxBoundary::new(100, 200, 0);
        let b = TxBoundary::new(150, 300, 0);
        let c = TxBoundary::new(200, 400, 0); // [200, 400) 不与 [100, 200) overlap（半开区间）
        let d = TxBoundary::new(50, 80, 0);

        assert!(a.overlaps(b)); // [100,200) ∩ [150,300) ✓
        assert!(b.overlaps(a)); // 对称
        assert!(!a.overlaps(c)); // [100,200) ∩ [200,400) ✗ (半开)
        assert!(!a.overlaps(d)); // [100,200) ∩ [50,80)   ✗
    }

    #[test]
    fn test_contains() {
        let outer = TxBoundary::new(100, 500, 0);
        let inner = TxBoundary::new(200, 400, 0);
        let partial = TxBoundary::new(200, 600, 0);

        assert!(outer.contains(inner));
        assert!(!outer.contains(partial));
        assert!(!inner.contains(outer));
    }

    #[test]
    fn test_sort_by_left() {
        let mut v = vec![
            TxBoundary::new(500, 600, 0),
            TxBoundary::new(100, 200, 0),
            TxBoundary::new(100, 150, 0),
        ];
        v.sort();
        assert_eq!(v[0].left(), 100);
        assert_eq!(v[0].right(), 150); // 同 left 时 right 小的在前
        assert_eq!(v[1].right(), 200);
        assert_eq!(v[2].left(), 500);
    }
}
