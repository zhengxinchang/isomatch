// use std::collections::HashSet;

use std::io::{self, Write};

use ahash::HashSet;

use crate::{
    classify::{
        class_code::*, query_ptir::QueryPTIR, ref_ptir::RefPTIR, ref_ptir_manager::RefPTIRManager,
    },
    core::{tx_strand::ISOMSTRAND, tx_type::TxType},
    index::fasta::FastaReader,
    merge::guide::GuideDb,
    utils::rev_comp,
};


type GeneId = u32;
#[derive(Debug, Clone)]
pub struct ClassifyRecord {
    // Group 0: query PTIR
    query_ptir: QueryPTIR,

    // Group 1: Query intrinsic properties (no external data needed)
    query_length: u32,
    query_exon_n: u16,
    query_strand: ISOMSTRAND,

    has_chr_in_ref_gtf: bool,

    per_junction_known_vec: Vec<bool>,
    all_junction_known: bool,

    per_junction_left_site_known: Vec<bool>,
    per_junction_right_site_known: Vec<bool>,
    all_splice_sites_known: bool,

    contained_in_known_intron: bool, //genic_intron

    has_intron_retention_against_catalog: bool, // overide the NIC subtyep with intron_retention

    same_strand_overlap_genes: Vec<GeneId>,

    antisense_genes: Vec<GeneId>,

    diff_to_gene_tss: i32,

    diff_to_gene_tts: i32,

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
    query_to_ref_matched_junctions: usize,
    query_to_ref_matched_exons: usize,
    bite: bool,

    // Group 3: Sequence-derived QC attributes (requires: reference genome FASTA)
    all_canonical: bool,
    _rts_stage: bool,
    perc_a_downstream_tts: Option<f32>,
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

impl ClassifyRecord {
    pub fn new(query_ptir: &QueryPTIR, ref_ptir_manager: &RefPTIRManager) -> Self {
        let mut cr = Self {
            // Group 0: query PTIR
            query_ptir: query_ptir.clone(),

            // Group 1: Query intrinsic properties (no external data needed)
            query_length: query_ptir.end() - query_ptir.start(),
            query_exon_n: query_ptir.n_exons(),
            query_strand: *query_ptir.standard(),

            has_chr_in_ref_gtf: true,

            per_junction_known_vec: Vec::new(),
            all_junction_known: false,

            per_junction_left_site_known: Vec::new(),
            per_junction_right_site_known: Vec::new(),
            all_splice_sites_known: false,

            contained_in_known_intron: false, //genic_intron

            has_intron_retention_against_catalog: false, // overide the NIC subtyep with intron_retention

            same_strand_overlap_genes: Vec::new(),

            antisense_genes: Vec::new(),

            diff_to_gene_tss: 0,

            diff_to_gene_tts: 0,

            // Group 2: Reference annotation relationship
            cc: ClassCode::Intergenic,
            ref_gene_id: String::new(),
            ref_gene_name: String::new(),
            ref_tx_id: String::new(),
            ref_length: 0,
            ref_exon_n: 0,
            ref_strand: ISOMSTRAND::Unknown,
            diff_to_tss: 0,
            diff_to_tes: 0,
            query_to_ref_matched_junctions: 0,
            query_to_ref_matched_exons: 0,
            bite: false,

            // Group 3: Sequence-derived QC attributes
            all_canonical: matches!(query_ptir.tx_type(), TxType::ALLC | TxType::MONO),
            _rts_stage: false,
            perc_a_downstream_tts: None,
            seq_a_downstream_tts: None,
            _poly_a_motif: None,
            _poly_a_motif_found: false,
            _poly_a_dist: None,

            // Group 4: Orthogonal evidence
            within_cage_peak: None,
            dist_to_cage_peak: None,
            within_poly_a_peak: None,
            dist_to_poly_a_site: None,
        };

        if query_ptir.n_exons() > 1 {
            if !ref_ptir_manager.has_chr(&query_ptir.chr_name) {
                cr.has_chr_in_ref_gtf = false;
                return cr;
            }

            let junctions = &query_ptir.junction_vec().as_ref().unwrap();
            {
                let (is_all_junction_match, junction_match_vec) = ref_ptir_manager.junction_match(
                    &query_ptir.chr_name,
                    &junctions,
                    query_ptir.standard(),
                );

                cr.has_chr_in_ref_gtf = true;
                cr.per_junction_known_vec = junction_match_vec;
                cr.all_junction_known = is_all_junction_match;
            }

            // update splice match
            {
                let (is_all_splice_sites_match, left_site_match_vec, right_site_match_vec) =
                    ref_ptir_manager.splice_site_match(&query_ptir.chr_name, &junctions);

                cr.has_chr_in_ref_gtf = true;
                cr.per_junction_left_site_known = left_site_match_vec;
                cr.per_junction_right_site_known = right_site_match_vec;
                cr.all_splice_sites_known = is_all_splice_sites_match;
            }

            // with in known junction
            cr.contained_in_known_intron = ref_ptir_manager.contained_in_known_intron(
                &query_ptir.chr_name,
                query_ptir.standard(),
                query_ptir.start(),
                query_ptir.end(),
            );

            // has_intron_retention_against_catalog
            cr.has_intron_retention_against_catalog =
                ref_ptir_manager.has_intron_retention_against_catalog(
                    &query_ptir.chr_name,
                    query_ptir.standard(),
                    &query_ptir.exons_vec(),
                );
        } else {
            todo!()
        }
        todo!()
    }


}
