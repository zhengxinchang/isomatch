//! struct of query PIIR

use crate::{
    core::{ptir::PTIR, string_pool::StringSpan, tx_base::TxBase},
    index::attributes_index::IsomSrcRecord,
    merge::policy::MergePolicyUsed,
};

pub struct QueryPTIR {
    pub is_isomatch_merged: bool,
    base: PTIR,
    src_vec: Vec<IsomSrcRecord>,
    merge_policies: (MergePolicyUsed, MergePolicyUsed, MergePolicyUsed),
}

impl QueryPTIR {}

pub struct QueryPTIRManager {}
