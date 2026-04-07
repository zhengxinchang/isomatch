use std::cmp::Ordering;

use crate::core::junction_pool::*;
use crate::core::splice_site_pool::{SpliceSitePool, SpliceSiteSpan};
use crate::core::string_pool::{StringPool, StringSpan};
use crate::core::tx_base_flag::TxBaseFlags;
use crate::core::{tx_base_error::TxBaseError, tx_boundary::TxBoundary};
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
    fn gtf_offset(&self) -> u64;
    fn gtf_len(&self) -> u32;
    fn n_exons(&self) -> u16;
    fn junctions(&self,junction_pool:&mut JunctionPool) -> Vec<(u32,u32)>;
    fn splice_sites(&self,splice_sites_pool:&mut SpliceSitePool) -> Vec<u8>;
    fn source_tx_id(&self,string_pool:&mut StringPool) -> String;
    fn source_gene_id(&self,string_pool:&mut StringPool) -> String;
    fn strand(&self) -> u8 {
        self.flags().get_strand()
    }
}

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
    pub _gtf_offset: u64, // byte offset of the GTF record in the original GTF file
    pub _gtf_len: u32,    // byte length of the GTF record in the original GTF file
    pub n_exons: u16,
    pub junctions: JunctionSpan,

    // new added
    pub splice_sites: SpliceSiteSpan,

    /// Direct reference into the on-disk string section for GTF `transcript_id`.
    pub tx_id_span: StringSpan,
    /// Direct reference into the on-disk string section for GTF `gene_id`.
    pub gene_id_span: StringSpan,
}

impl TxBase {
    pub fn new(
        tx_idx: u32, // record index in the GTF file
        chrom_id: u16,
        start: u32,
        end: u32,
        strand: u8,
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
            tx_idx: tx_idx,
            boundary: TxBoundary::new(start, end, strand),
            chrom_id,
            start,
            end,
            flags: TxBaseFlags::new(strand, seq_hash != 0)?,
            seq_hash,
            ref_hash,
            _gtf_offset: 0,
            _gtf_len: 0,
            n_exons,
            splice_sites: splice_site_span,
            junctions: junction_span,
            tx_id_span: transcript_span,
            gene_id_span: gene_span,
        })
    }

    pub fn strand(&self) -> u8 {
        self.flags.get_strand()
    }

    pub fn sort_key(&self) -> (u16, u32, u32, u8) {
        (self.chrom_id, self.start, self.end, self.strand())
    }

    pub fn junction_slice<'a>(&self, pool: &'a JunctionPool) -> Result<&'a [u32], TxBaseError> {
        pool.get(self.junctions)
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
