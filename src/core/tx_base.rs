use std::cmp::Ordering;

use crate::core::junction_pool::*;
use crate::core::splice_site_pair::SpliceSitePair;
use crate::core::splice_site_pool::SpliceSitePool;
use crate::core::splice_site_span::SpliceSiteSpan;
use crate::core::string_pool::{StringPool, StringSpan};
use crate::core::tx_base_flag::TxBaseFlags;
use crate::core::tx_strand::ISOMSTRAND;
use crate::core::{core_error::TxBaseError, tx_boundary::TxBoundary};
pub trait TxBaseTrait {
    fn tx_idx(&self) -> u32;
    fn tx_boundary(&self) -> TxBoundary {
        TxBoundary::new(self.start(), self.end(), self.strand())
    }
    fn chrom_id(&self) -> u16;
    fn start(&self) -> u32;
    fn end(&self) -> u32;
    fn flags(&self) -> TxBaseFlags;
    fn seq_hash(&self) -> u128;
    fn ref_hash(&self) -> u128;
    // fn gtf_offset(&self) -> u64;
    // fn gtf_len(&self) -> u32;
    fn n_exons(&self) -> u16;
    fn junctions(&self, junction_pool: &JunctionPool, string_pool: &StringPool) -> Vec<(u32, u32)>;
    fn splice_sites(
        &self,
        splice_sites_pool: &SpliceSitePool,
        string_pool: &StringPool,
    ) -> Vec<SpliceSitePair>;
    fn source_tx_id(&self, string_pool: &StringPool) -> String;
    fn source_gene_id(&self, string_pool: &StringPool) -> String;
    fn strand(&self) -> ISOMSTRAND {
        self.flags().get_strand()
    }
}

/// core data stucture for transcript
/// for persistance on disk.
#[derive(Debug, Clone, Copy, Hash)]
pub struct TxBase {
    pub tx_idx: u32,
    pub boundary: TxBoundary,
    pub chrom_id: u16,
    pub start: u32,
    pub end: u32,
    pub flags: TxBaseFlags,
    pub seq_hash: u128,
    pub ref_hash: u128,
    pub n_exons: u16,
    pub junctions_span: JunctionSpan,

    // new added
    pub splice_sites_span: SpliceSiteSpan,

    /// Direct reference into the on-disk string section for GTF `transcript_id`.
    pub tx_id_span: StringSpan,
    /// Direct reference into the on-disk string section for GTF `gene_id`.
    pub gene_id_span: StringSpan,
}

impl TxBase {
    pub fn new(
        gid: u32,
        chrom_id: u16,
        start: u32,
        end: u32,
        strand: ISOMSTRAND,
        seq_hash: u128,
        ref_hash: u128,
        n_exons: u16,
        splice_site_span: SpliceSiteSpan,
        junction_span: JunctionSpan,
        transcript_span: StringSpan,
        gene_span: StringSpan,
    ) -> Result<Self, TxBaseError> {
        if start > end {
            return Err(TxBaseError::InvalidBounds { start, end });
        }
        if n_exons == 0 {
            return Err(TxBaseError::InvalidExonCount { n_exons });
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

    pub fn strand(&self) -> ISOMSTRAND {
        self.flags.get_strand()
    }

    pub fn sort_key(&self) -> (u16, u32, u32, ISOMSTRAND) {
        (self.chrom_id, self.start, self.end, self.strand())
    }

    pub fn junction_slice<'a>(&self, pool: &'a JunctionPool) -> Result<&'a [u32], TxBaseError> {
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
