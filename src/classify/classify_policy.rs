use crate::classify::{classification_code::*, query_ptir::QueryPTIR, ref_ptir::RefPTIRManager};

pub struct Classification {
    cc: ClassCode,
    dis_to_tss: i32,
    dis_to_tes: i32,
    ref_gene_id: String,
    ref_gene_name: String,
    ref_tx_id: String,
}

pub fn classify(query_ptir: &QueryPTIR, ref_ptir: &RefPTIRManager) -> Classification {
    todo!()
}
