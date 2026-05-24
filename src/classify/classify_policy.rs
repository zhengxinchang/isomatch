use crate::{
    classify::{
        class_code::*,
        query_ptir::QueryPTIR,
        ref_ptir::{RefPTIR, RefPTIRManager},
    },
    core::tx_strand::ISOMSTRAND,
    index::fasta::FastaReader,
    merge::guide::GuideDb,
};

pub struct Classification {
    // Group 1: Query intrinsic properties (no external data needed)
    query_length: u32,
    query_exon_n: u32,
    query_strand: ISOMSTRAND,

    // Group 2: Reference annotation relationship (requires: query GTF + reference GTF)
    cc: ClassCode,
    ref_gene_id: String,
    ref_gene_name: String,
    ref_tx_id: String,
    ref_length: u32,
    ref_exon_n: u32,
    ref_strand: ISOMSTRAND,
    diff_to_tss: i32,
    diff_to_tes: i32,
    _diff_to_gene_tss: i32,
    _diff_to_gene_tts: i32,
    query_to_ref_matched_junctions: usize,
    query_to_ref_matched_exons: usize,
    bite: bool,

    // Group 3: Sequence-derived QC attributes (requires: reference genome FASTA)
    all_canonical: bool,
    _rts_stage: bool,
    perc_a_downstream_tts: f32,
    seq_a_downstream_tts: Option<String>,
    _poly_a_motif: Option<String>,
    _poly_a_motif_found: bool,
    _poly_a_dist: Option<i32>,

    // Group 4: Orthogonal evidence (requires: external BED files, optional)
    within_cage_peak: Option<bool>,
    dist_to_cage_peak: Option<i32>,
    within_poly_a_peak: Option<bool>,
    dist_to_poly_a_site: Option<i32>,
}

impl Classification {
    pub fn new(query_ptir: &QueryPTIR) -> Self {
        // create new struct and update group 1
        todo!()
    }
}

pub fn classify(
    query_ptir: &QueryPTIR,
    ref_ptir: &RefPTIR,
    ref_tss: Option<&GuideDb>,
    ref_tes: Option<&GuideDb>,
    ref_fa: &FastaReader,
) -> Classification {
    let mut class = Classification::new(&query_ptir);
    let cc = get_class_code(query_ptir, ref_ptir);
    class.cc = cc;

    // update group2 based on query and refernece ptir

    // update group3 and group4 based on /ssd2/projects/isomatch-dev/src/classify/group3and4.txt

    // return Classification
    todo!()
}

pub fn get_class_code(query_ptir: &QueryPTIR, ref_ptir: &RefPTIR) -> ClassCode {
    // get a classificiaton code base on query and ref
    // follow up squanti3 rule in /ssd2/projects/isomatch-dev/src/classify/class_code.rs
    todo!()
}
