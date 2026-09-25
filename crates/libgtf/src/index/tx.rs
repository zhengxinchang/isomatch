use std::cmp::Ordering;

use crate::error::IndexDataError;
use crate::gtf::Strand;
use crate::index::junction_pool::{JunctionPool, JunctionSpan};
use crate::index::splice_site_pool::SpliceSitePair;
use crate::index::splice_site_pool::SpliceSitePool;
use crate::index::splice_site_pool::SpliceSiteSpan;
use crate::index::string_pool::{StringPool, StringSpan};
use crate::traits::{DiskSize, Encodable, PartialLoad};

// pub trait TxBaseTrait {
//     fn tx_idx(&self) -> u64;
//     fn tx_boundary(&self) -> TxBoundary {
//         TxBoundary::new(self.start(), self.end(), self.strand())
//     }
//     fn chrom_id(&self) -> u16;
//     fn start(&self) -> u32;
//     fn end(&self) -> u32;
//     fn flags(&self) -> TxBaseFlags;
//     fn seq_hash(&self) -> u128;
//     fn ref_hash(&self) -> u128;
//     // fn gtf_offset(&self) -> u64;
//     // fn gtf_len(&self) -> u32;
//     fn n_exons(&self) -> u16;
//     fn junctions(&self, junction_pool: &JunctionPool, string_pool: &StringPool) -> Vec<(u32, u32)>;
//     fn splice_sites(
//         &self,
//         splice_sites_pool: &SpliceSitePool,
//         string_pool: &StringPool,
//     ) -> Vec<SpliceSitePair>;
//     fn source_tx_id(&self, string_pool: &StringPool) -> String;
//     fn source_gene_id(&self, string_pool: &StringPool) -> String;
//     fn strand(&self) -> Strand {
//         self.flags().get_strand()
//     }
// }

/// core data stucture for transcript
/// for persistance on disk.
#[derive(Debug, Clone, Copy)]
pub struct TxBase {
    tx_idx: u64,
    boundary: TxBoundary,
    chrom_id: u16,
    start: u32,
    end: u32,
    flags: TxBaseFlags,
    seq_hash: u128,
    ref_hash: u128,
    n_exons: u16,
    junctions_span: JunctionSpan,

    // new added
    splice_sites_span: SpliceSiteSpan,

    /// Direct reference into the on-disk string section for GTF `transcript_id`.
    tx_id_span: StringSpan,
    /// Direct reference into the on-disk string section for GTF `gene_id`.
    gene_id_span: StringSpan,
}

impl TxBase {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        gid: u64,
        chrom_id: u16,
        start: u32,
        end: u32,
        strand: Strand,
        seq_hash: u128,
        ref_hash: u128,
        n_exons: u16,
        splice_site_span: SpliceSiteSpan,
        junction_span: JunctionSpan,
        transcript_span: StringSpan,
        gene_span: StringSpan,
    ) -> Result<Self, IndexDataError> {
        if start > end {
            return Err(IndexDataError::InvalidBounds { start, end });
        }
        if n_exons == 0 {
            return Err(IndexDataError::InvalidExonCount { n_exons });
        }

        Ok(Self {
            tx_idx: gid,
            boundary: TxBoundary::new(start, end, strand),
            chrom_id,
            start,
            end,
            flags: TxBaseFlags::new(strand, seq_hash != 0)?,
            seq_hash,
            ref_hash,
            n_exons,
            splice_sites_span: splice_site_span,
            junctions_span: junction_span,
            tx_id_span: transcript_span,
            gene_id_span: gene_span,
        })
    }

    pub fn strand(&self) -> Strand {
        self.flags.get_strand()
    }

    pub fn sort_key(&self) -> (u16, u32, u32, Strand) {
        (self.chrom_id, self.start, self.end, self.strand())
    }

    pub fn junction_slice<'a>(&self, pool: &'a JunctionPool) -> Result<&'a [u32], IndexDataError> {
        pool.get(self.junctions_span)
    }
}

impl PartialOrd for TxBase {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for TxBase {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.start, self.end, self.strand()).cmp(&(other.start, other.end, other.strand()))
    }
}

impl Eq for TxBase {}

impl PartialEq for TxBase {
    fn eq(&self, other: &Self) -> bool {
        self.start == other.start && self.end == other.end && (self.strand() == other.strand())
    }
}

/// Flags for TxBase.
/// bit 0-1: strand (0 for +, 1 for -, 2 for unknown)
/// bit 2: seq_hash is valid
/// bit 3: is gtf output of isomatch merge? 1 for yes, 0 for no.
/// bit 3bit * 4 = 12 bit.  merge policy, SJ TSS TES MONO, each is guided, major, longer, shorter, [place holder] * 4
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
#[repr(transparent)]
pub struct TxBaseFlags(pub u16);

impl TxBaseFlags {
    const STRAND_MASK: u16 = 0b11;
    const HAS_SEQ_HASH_BIT: u16 = 1 << 2;

    pub fn new(strand: Strand, seq_has_hash: bool) -> Result<Self, IndexDataError> {
        let mut flags = Self(u16::from(strand.to_bit()));

        if seq_has_hash {
            flags.0 |= Self::HAS_SEQ_HASH_BIT;
        }

        Ok(flags)
    }

    pub fn get_strand(self) -> Strand {
        Strand::from_bit((self.0 & Self::STRAND_MASK) as u8).unwrap()
    }

    pub fn set_strand(&mut self, strand: Strand) -> Result<(), IndexDataError> {
        self.0 &= !Self::STRAND_MASK;
        self.0 |= u16::from(strand.to_bit());
        Ok(())
    }

    pub fn get_seq_has_hash(&self) -> bool {
        self.0 & Self::HAS_SEQ_HASH_BIT == Self::HAS_SEQ_HASH_BIT
    }

    pub fn set_seq_has_hash(&mut self, has_hash: bool) {
        if has_hash {
            self.0 |= Self::HAS_SEQ_HASH_BIT;
        } else {
            self.0 &= !Self::HAS_SEQ_HASH_BIT;
        }
    }

    pub fn bits(self) -> u16 {
        self.0
    }

    // pub fn strand(self) -> u8 {
    //     self.get_strand()
    // }

    pub fn seq_has_hash(self) -> bool {
        self.get_seq_has_hash()
    }
}

/// Encodes a transcript boundary as a single u64:
/// [left: 32bit | right: 30bit | strand: 2bit]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(transparent)]
pub struct TxBoundary(u64);

impl TxBoundary {
    const STRAND_MASK: u64 = 0b11;
    const RIGHT_MASK: u64 = 0x3FFF_FFFF;

    #[inline(always)]
    pub fn new(left: u32, right: u32, strand: Strand) -> Self {
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
    pub fn strand(self) -> Strand {
        Strand::from_bit((self.0 & Self::STRAND_MASK) as u8).unwrap()
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
            Strand::Plus => '+',
            Strand::Minus => '-',
            Strand::Unknown => '.',
        };
        write!(f, "[{}, {}]{}", self.left(), self.right(), strand_ch)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct TxBaseLoadArgs {
    pub chrom_id: u16,
}

impl TxBase {
    // 标注一下这个不应该被使用
    pub fn tx_idx(&self) -> u64 {
        self.tx_idx
    }
    pub fn tx_boundary(&self) -> TxBoundary {
        self.boundary
    }
    pub fn chrom_id(&self) -> u16 {
        self.chrom_id
    }

    pub fn start(&self) -> u32 {
        self.start
    }

    pub fn end(&self) -> u32 {
        self.end
    }

    pub fn flags(&self) -> TxBaseFlags {
        self.flags
    }

    pub fn seq_hash(&self) -> u128 {
        self.seq_hash
    }

    pub fn ref_hash(&self) -> u128 {
        self.ref_hash
    }

    // fn gtf_offset(&self) -> u64 {
    //     0
    // }

    // fn gtf_len(&self) -> u32 {
    //     0
    // }

    pub fn n_exons(&self) -> u16 {
        self.n_exons
    }

    pub fn junctions(
        &self,
        junction_pool: &JunctionPool,
        string_pool: &StringPool,
    ) -> Vec<(u32, u32)> {
        let raw = junction_pool.get(self.junctions_span).unwrap();

        let chunks = raw.chunks_exact(2);
        let remainder = chunks.remainder();
        assert!(
            remainder.is_empty(),
            "junction coordinate count must be even for transcript {}, got {}",
            self.source_tx_id(string_pool),
            raw.len()
        );

        chunks.map(|pair| (pair[0], pair[1])).collect()
    }

    pub fn splice_sites(&self, splice_sites_pool: &SpliceSitePool) -> Vec<SpliceSitePair> {
        splice_sites_pool
            .get_pair(self.splice_sites_span)
            .unwrap()
            .to_vec()
    }

    pub fn source_tx_id(&self, string_pool: &StringPool) -> String {
        string_pool.get(self.tx_id_span).unwrap().to_owned()
    }

    pub fn source_gene_id(&self, string_pool: &StringPool) -> String {
        string_pool.get(self.gene_id_span).unwrap().to_owned()
    }
}

impl DiskSize for TxBase {
    // Current layout drops on-disk chrom_id and the unused gtf_offset/gtf_len fields.
    const DISK_SIZE: usize = 96;
}

impl Encodable for TxBase {
    type Error = IndexDataError;

    /// no need to do the TxBoundary encoding here since TxBase already stores start, end, strand separately for easy access
    fn encode_to<W: std::io::Write>(&self, writer: &mut W) -> Result<usize, Self::Error> {
        writer
            .write_all(&self.tx_idx.to_le_bytes())
            .map_err(|e| IndexDataError::Io(e.to_string()))?;
        writer
            .write_all(&self.start.to_le_bytes())
            .map_err(|e| IndexDataError::Io(e.to_string()))?;
        writer
            .write_all(&self.end.to_le_bytes())
            .map_err(|e| IndexDataError::Io(e.to_string()))?;
        writer
            .write_all(&self.flags.bits().to_le_bytes())
            .map_err(|e| IndexDataError::Io(e.to_string()))?;
        writer
            .write_all(&self.seq_hash.to_le_bytes())
            .map_err(|e| IndexDataError::Io(e.to_string()))?;
        writer
            .write_all(&self.ref_hash.to_le_bytes())
            .map_err(|e| IndexDataError::Io(e.to_string()))?;
        writer
            .write_all(&self.n_exons.to_le_bytes())
            .map_err(|e| IndexDataError::Io(e.to_string()))?;
        writer
            .write_all(&self.junctions_span.offset.to_le_bytes())
            .map_err(|e| IndexDataError::Io(e.to_string()))?;
        writer
            .write_all(&self.junctions_span.count.to_le_bytes())
            .map_err(|e| IndexDataError::Io(e.to_string()))?;
        writer
            .write_all(&self.splice_sites_span.offset.to_le_bytes())
            .map_err(|e| IndexDataError::Io(e.to_string()))?;
        writer
            .write_all(&self.splice_sites_span.count.to_le_bytes())
            .map_err(|e| IndexDataError::Io(e.to_string()))?;
        writer
            .write_all(&self.tx_id_span.offset.to_le_bytes())
            .map_err(|e| IndexDataError::Io(e.to_string()))?;
        writer
            .write_all(&self.tx_id_span.byte_len.to_le_bytes())
            .map_err(|e| IndexDataError::Io(e.to_string()))?;
        writer
            .write_all(&self.gene_id_span.offset.to_le_bytes())
            .map_err(|e| IndexDataError::Io(e.to_string()))?;
        writer
            .write_all(&self.gene_id_span.byte_len.to_le_bytes())
            .map_err(|e| IndexDataError::Io(e.to_string()))?;
        Ok(Self::DISK_SIZE)
    }
}

impl PartialLoad for TxBase {
    type Error = IndexDataError;
    type Args = TxBaseLoadArgs;

    fn load_range<R: std::io::Read + std::io::Seek>(
        reader: &mut R,
        offset: u64,
        _len: usize, // always DISK_SIZE for fixed-size TxBase, ignored
        args: Self::Args,
    ) -> Result<Self, Self::Error> {
        reader
            .seek(std::io::SeekFrom::Start(offset))
            .map_err(|e| IndexDataError::Io(e.to_string()))?;

        let mut buf = [0u8; TxBase::DISK_SIZE];
        reader
            .read_exact(&mut buf)
            .map_err(|e| IndexDataError::Io(e.to_string()))?;

        let tx_id = u64::from_le_bytes(buf[0..8].try_into().unwrap());
        let start = u32::from_le_bytes(buf[8..12].try_into().unwrap());
        let end = u32::from_le_bytes(buf[12..16].try_into().unwrap());
        let flags = TxBaseFlags(u16::from_le_bytes(buf[16..18].try_into().unwrap()));
        let seq_hash = u128::from_le_bytes(buf[18..34].try_into().unwrap());
        let ref_hash = u128::from_le_bytes(buf[34..50].try_into().unwrap());
        let n_exons = u16::from_le_bytes(buf[50..52].try_into().unwrap());
        let junctions_offset = u32::from_le_bytes(buf[52..56].try_into().unwrap());
        let junctions_count = u16::from_le_bytes(buf[56..58].try_into().unwrap());
        let splice_sites_offset = u32::from_le_bytes(buf[58..62].try_into().unwrap());
        let splice_sites_count = u16::from_le_bytes(buf[62..64].try_into().unwrap());
        let transcript_span_offset = u64::from_le_bytes(buf[64..72].try_into().unwrap());
        let transcript_span_byte_len = u64::from_le_bytes(buf[72..80].try_into().unwrap());
        let gene_span_offset = u64::from_le_bytes(buf[80..88].try_into().unwrap());
        let gene_span_byte_len = u64::from_le_bytes(buf[88..96].try_into().unwrap());

        Ok(Self {
            tx_idx: tx_id,
            boundary: TxBoundary::new(start, end, flags.get_strand()),
            chrom_id: args.chrom_id,
            start,
            end,
            flags,
            seq_hash,
            ref_hash,
            n_exons,
            junctions_span: JunctionSpan {
                offset: junctions_offset,
                count: junctions_count,
            },
            splice_sites_span: SpliceSiteSpan {
                offset: splice_sites_offset,
                count: splice_sites_count,
            },
            tx_id_span: StringSpan {
                offset: transcript_span_offset,
                byte_len: transcript_span_byte_len,
            },
            gene_id_span: StringSpan {
                offset: gene_span_offset,
                byte_len: gene_span_byte_len,
            },
        })
    }
}
