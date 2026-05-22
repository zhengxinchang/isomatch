//! struct of query PIIR

use crate::{
    classify::isofrom_src_entry::IsoSrcSpan,
    core::{ptir::PTIR, tx_base::TxBase},
};

pub struct QueryPTIR {
    txbase: PTIR,
    isoform_src_span: IsoSrcSpan,
}
