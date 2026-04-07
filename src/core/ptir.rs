//! PTIR refers to Pan-transcript intermediate representation
//! This object is the core data structure that is loaded from
//! index file.
//! it is used in the merge and annotate command.

use crate::core::{
    junction_pool::JunctionPool,
    splice_site_pool::{SpliceSitePair, SpliceSitePool},
    status::TxType,
    string_pool::StringPool,
    tx_base::{TxBase, TxBaseTrait},
    tx_boundary::TxBoundary,
};

pub struct PTIR {
    pub tx_boundary: TxBoundary,
    pub start: u32,
    pub end: u32,
    pub strand: u8,
    pub n_exons: u16,
    pub refhash: u128,
    pub seqhash: Option<u128>,
    pub junction_vec: Option<Vec<(u32, u32)>>,
    pub splice_site_vec: Option<Vec<SpliceSitePair>>,
    pub tx_type: TxType,
    pub source_txid: String,
    pub source_geneid: String,
}

impl PTIR {
    pub fn from_tx_base(
        tb: TxBase,
        junc_pool: &JunctionPool,
        spl_site_pool: &SpliceSitePool,
        string_pool: &StringPool,
    ) -> Self {
        let splice_site_vec: Option<Vec<SpliceSitePair>> = if tb.n_exons() == 1 {
            Some(tb.splice_sites(spl_site_pool))
        } else {
            None
        };

        let tx_type = match &splice_site_vec {
            Some(vec) => {
                if vec.iter().all(|&b| b.is_canonical()) {
                    TxType::ALLC
                } else if !vec.iter().all(|&b| b.is_canonical()) {
                    TxType::NOTC
                } else {
                    TxType::PRTC
                }
            }
            None => TxType::MONO,
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
            junction_vec: if tb.n_exons() == 1 {
                Some(tb.junctions(junc_pool))
            } else {
                None
            },
            splice_site_vec: splice_site_vec,
            tx_type: tx_type,
            source_txid: tb.source_gene_id(string_pool),
            source_geneid: tb.source_tx_id(string_pool),
        }
    }
}


impl PartialEq for PTIR {
    fn eq(&self, other: &Self) -> bool {
        self.tx_boundary == other.tx_boundary 
    }
}

impl Eq for PTIR {

}

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