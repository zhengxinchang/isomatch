//! PTIR refers to Pan-transcript intermediate representation
//! This object is the core data structure that is loaded from
//! index file.
//! it is used in the merge and annotate command.

use crate::core::{
    junction_pool::JunctionPool,
    splice_site_pair::SpliceSitePair,
    splice_site_pool::SpliceSitePool,
    string_pool::StringPool,
    tx_base::{TxBase, TxBaseTrait},
    tx_boundary::TxBoundary,
    tx_strand::ISOMSTRAND,
    tx_type::TxType,
};

pub struct PTIR {
    pub tx_boundary: TxBoundary,
    pub start: u32,
    pub end: u32,
    pub strand: ISOMSTRAND,
    pub n_exons: u16,
    pub refhash: u128,
    pub seqhash: Option<u128>,
    pub junction_vec: Option<Vec<(u32, u32)>>,
    pub splice_site_vec: Option<Vec<SpliceSitePair>>,
    pub tx_type: TxType,
    pub source_file_id: usize,
    pub source_txid: String,
    pub source_geneid: String,
}

impl PTIR {
    pub fn from_tx_base(
        tb: TxBase,
        file_id: usize,
        junc_pool: &JunctionPool,
        spl_site_pool: &SpliceSitePool,
        string_pool: &StringPool,
    ) -> Self {
        let splice_site_vec: Option<Vec<SpliceSitePair>> = if tb.n_exons() > 1 {
            Some(tb.splice_sites(spl_site_pool))
        } else {
            None
        };

        let tx_type = match splice_site_vec.as_ref() {
            None => TxType::MONO,
            Some(vec) if vec.iter().all(|site| site.is_canonical()) => TxType::ALLC,
            Some(vec) if vec.iter().any(|site| site.is_canonical()) => TxType::PRTC,
            Some(_) => TxType::NOTC,
        };

        Self {
            tx_boundary: tb.tx_boundary(),
            start: tb.start(),
            end: tb.end(),
            strand: tb.strand(),
            n_exons: tb.n_exons(),
            refhash: tb.ref_hash(),
            seqhash: if tb.flags.get_seq_has_hash() {
                Some(tb.seq_hash())
            } else {
                None
            },
            junction_vec: if tb.n_exons() > 1 {
                Some(tb.junctions(junc_pool))
            } else {
                None
            },
            splice_site_vec: splice_site_vec,
            tx_type: tx_type,
            source_file_id: file_id,
            source_txid: tb.source_tx_id(string_pool),
            source_geneid: tb.source_gene_id(string_pool),
        }
    }

    pub fn overlap(&self, other: &PTIR) -> bool {
        self.tx_boundary.overlaps(other.tx_boundary)
    }

    pub fn junctions(&self) -> Option<&[(u32, u32)]> {
        self.junction_vec.as_deref()
    }

    pub fn tss(&self) -> u32 {
        match self.strand {
            ISOMSTRAND::Plus => self.start,
            ISOMSTRAND::Minus => self.end,
            ISOMSTRAND::Unknown => self.start,
        }
    }

    pub fn tes(&self) -> u32 {
        match self.strand {
            ISOMSTRAND::Plus => self.end,
            ISOMSTRAND::Minus => self.start,
            ISOMSTRAND::Unknown => self.end,
        }
    }
}

impl PartialEq for PTIR {
    fn eq(&self, other: &Self) -> bool {
        self.tx_boundary == other.tx_boundary
    }
}

impl Eq for PTIR {}

impl PartialOrd for PTIR {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PTIR {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.tx_boundary.cmp(&other.tx_boundary)
    }
}

impl std::fmt::Display for PTIR {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}-{}, strand={}, n_exons={}, tx_id={}, gene_id={}, file_id={}",
            self.start,
            self.end,
            self.strand,
            self.n_exons,
            self.source_txid,
            self.source_geneid,
            self.source_file_id
        )
    }
}
