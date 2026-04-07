//! PTIR refers to Pan-transcript intermediate representation
//! This object is the core data structure that is loaded from 
//! index file. 
//! it is used in the merge and annotate command.
//! 
//! 

use crate::core::{status::TxType, tx_base::{TxBase, TxBaseTrait}, tx_boundary::TxBoundary};

pub struct PTIR {
    pub global_chrom_id:u32,
    pub tx_boundary: TxBoundary,
    pub start:u32,
    pub end:u32,
    pub strand:u8,
    pub n_exons:u16,
    pub refhash:u128,
    pub seqhash:Option<u128>,
    pub junction_vec: Option<Vec<(u32,u32)>>,
    pub splice_site_vec:Option<Vec<u8>>,
    pub tx_type:TxType,
    pub source_txid:String,
    pub source_geneid:String,
}

impl PTIR {

    pub fn from_tx_base(tx_base:TxBase) -> Self {
        let _ = tx_base;
        unimplemented!("PTIR::from_tx_base is not implemented yet")
    }

}
