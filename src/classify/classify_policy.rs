use std::io::{self, Write};

use ahash::{HashMap, HashSet};
use log::warn;

use crate::{
    ClassifyArgs,
    classify::{
        class_code::*,
        compare::{
            JunctionMatch, classify_junction_chain, exon_overlap_bases,
            find_consecutive_junction_chain, genomic_overlap_bases,
            mono_query_contained_in_ref_exon, mono_query_spans_ref_intron,
            query_has_intron_retention_against_ref, same_strand_transcript_space_tss_tes_diffs,
            splice_site_agreement,
        },
        query_ptir::QueryPTIR,
        ref_ptir::RefPTIR,
        ref_ptir_manager::RefPTIRManager,
    },
    constants::MOTIFS,
    core::{tx_strand::ISOMSTRAND, tx_type::TxType},
    index::fasta::FastaReader,
    merge::guide::GuideDb,
    utils::rev_comp,
};

/// Early per-reference classification before gene/catalog refinement.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PreClass {
    /// Candidate is already a full-splice match.
    Fsm(SubFSM),
    /// Candidate is already an incomplete-splice match.
    Ism(SubISM),
    // /// Candidate is already a direct novel-in-catalog mono-exon/intron-retention case.
    Nic(SubNIC),
    /// Candidate shares at least one known junction pair but is not FSM/ISM.
    AnyKnownJunction,
    /// Candidate shares at least one splice site but no full junction pair.
    AnyKnownSpliceSite,
    /// Candidate overlaps exonic sequence but no junction/splice-site evidence matched.
    GeneOverlap,
    /// Candidate has no useful same-strand structural relationship.
    None,
}

impl PreClass {
    /// Ranking used to keep the best reference transcript within each gene.
    fn rank(self) -> u8 {
        match self {
            Self::Fsm(_) => 6,
            Self::Ism(_) => 5,
            Self::Nic(_) => 4,
            Self::AnyKnownJunction => 3,
            Self::AnyKnownSpliceSite => 2,
            Self::GeneOverlap => 1,
            Self::None => 0,
        }
    }
}

/// One query-vs-reference transcript hit with the measurements needed for ranking.
/// One query may have multiple of this.
#[derive(Debug, Clone)]
struct CandidateHit {
    query_ptir: QueryPTIR,
    ref_ptir: RefPTIR,
    pre_class: PreClass,
    gene_id: String,
    gene_name: String,
    tx_id: String,
    ref_start: u32,
    ref_end: u32,
    ref_length: u32,
    ref_exons: u16,
    ref_strand: ISOMSTRAND,
    diff_tss: Option<i32>,
    diff_tes: Option<i32>,
    splice_site_hits: usize,
    matched_junctions: usize,
    exon_overlap_bases: u32,
    genomic_overlap_bases: u32,
}

impl CandidateHit {
    /// Build a candidate hit from one reference transcript plus an already-computed pre-class.
    ///
    fn from_ref(query: &QueryPTIR, reference: &RefPTIR, pre_class: PreClass) -> Self {
        // Cache exon lists once because several evidence metrics reuse them.
        let query_exons = query.exons_vec();
        let ref_exons = reference.exons_vec();
        let (diff_tss, diff_tes) = same_strand_transcript_space_tss_tes_diffs(query, reference);

        // get number of matched junctions.
        let matched_junctions = match (
            query.junction_vec().as_deref(),
            reference.junction_vec().as_deref(),
        ) {
            (Some(query_junctions), Some(ref_junctions)) => {
                let ref_set: HashSet<(u32, u32)> = ref_junctions.iter().copied().collect();
                query_junctions
                    .iter()
                    .filter(|junction| ref_set.contains(junction))
                    .count()
            }
            _ => 0,
        };

        Self {
            query_ptir: query.clone(),
            ref_ptir: reference.clone(),
            pre_class,
            gene_id: reference.base.source_geneid.clone(),
            gene_name: reference.reference_gene_name().to_string(),
            tx_id: match pre_class {
                PreClass::Fsm(_) | PreClass::Ism(_) => reference.base.source_txid.clone(),
                _ => "novel".to_string(),
            },
            ref_start: reference.start(),
            ref_end: reference.end(),
            ref_length: reference.transcript_len(),
            ref_exons: reference.n_exons(),
            ref_strand: *reference.standard(),
            diff_tss: matches!(pre_class, PreClass::Fsm(_) | PreClass::Ism(_)).then_some(diff_tss),
            diff_tes: matches!(pre_class, PreClass::Fsm(_) | PreClass::Ism(_)).then_some(diff_tes),
            splice_site_hits: splice_site_agreement(&query_exons, &ref_exons),
            matched_junctions,
            exon_overlap_bases: exon_overlap_bases(&query_exons, &ref_exons),
            genomic_overlap_bases: genomic_overlap_bases(
                query.start(),
                query.end(),
                reference.start(),
                reference.end(),
            ),
        }
    }

    /// Return true when `self` should replace another hit for the same gene.

    fn is_better_than(&self, other: &Self, query_exon_n: u16, query_len: u32) -> bool {
        let self_rank = self.pre_class.rank();
        let other_rank = other.pre_class.rank();
        if self_rank != other_rank {
            return self_rank > other_rank;
        }

        // must consider the other's preclass

        match self.pre_class {
            PreClass::Fsm(_) | PreClass::Ism(_) => {
                // make sure other pre_class is also FSM/ISM for comparison with endopoint difference tie-breaking.
                // otherwise, self always wins because it is FSM/ISM
                match other.pre_class {
                    PreClass::Fsm(_) | PreClass::Ism(_) => {
                        self.endpoint_total_diff() < other.endpoint_total_diff()
                    }
                    _ => return true,
                }
            }
            PreClass::AnyKnownJunction => {
                self.splice_site_hits > other.splice_site_hits
                    || (self.splice_site_hits == other.splice_site_hits
                        && self.exon_overlap_bases > other.exon_overlap_bases)
                    || (self.splice_site_hits == other.splice_site_hits
                        && self.exon_overlap_bases == other.exon_overlap_bases
                        && self.exon_count_diff(query_exon_n) < other.exon_count_diff(query_exon_n))
            }
            _ => self.score(query_exon_n, query_len) > other.score(query_exon_n, query_len),
        }
    }

    /// Total absolute endpoint difference for FSM/ISM tie-breaking.
    fn endpoint_total_diff(&self) -> i32 {
        // if !matches!(self.pre_class, PreClass::Fsm(_) | PreClass::Ism(_)) {
        //     warn!(
        //         "endpoint_total_diff called on non-FSM/ISM hit: treating as 0 for tie-breaking. preclass {:?} Query: {}, Ref: {}, Query PTIR: {:#?}, Ref PTIR: {:#?}",
        //         &self.pre_class, self.tx_id, self.gene_id, self.query_ptir, self.ref_ptir
        //     );
        //     return 0;
        // }

        // let diff_tss = self.diff_tss.unwrap_or(0);
        // let diff_tes = self.diff_tes.unwrap_or(0);
        // diff_tss.abs() + diff_tes.abs()

        match (self.diff_tss, self.diff_tes) {
            (Some(diff_tss), Some(diff_tes)) => diff_tss.abs() + diff_tes.abs(),
            _ => i32::MAX,
        }
    }

    /// Difference in exon count between query and this reference candidate.
    fn exon_count_diff(&self, query_exon_n: u16) -> u16 {
        query_exon_n.abs_diff(self.ref_exons)
    }

    /// Composite ranking score used after class-rank ties are handled.

    fn score(&self, query_exon_n: u16, query_len: u32) -> f64 {
        let query_len = query_len.max(1) as f64;
        let ref_len = self.ref_length.max(1) as f64;
        self.pre_class.rank() as f64
            + self.splice_site_hits as f64
            + self.exon_overlap_bases as f64 / query_len
            + self.genomic_overlap_bases as f64 / ref_len
            - self.exon_count_diff(query_exon_n) as f64
    }
}

/// Final per-query classification record written to the classification table.
#[derive(Debug, Clone)]
pub struct ClassifyRecord {
    query_ptir: QueryPTIR,

    query_length: u32,
    query_exon_n: u16,
    query_strand: ISOMSTRAND,

    has_chr_in_ref_gtf: bool,
    per_junction_known_vec: Vec<bool>,
    all_junction_known: bool,
    per_junction_left_site_known: Vec<bool>,
    per_junction_right_site_known: Vec<bool>,
    all_splice_sites_known: bool,
    contained_in_known_intron: bool,
    has_intron_retention_against_catalog: bool,
    same_strand_overlap_genes: Vec<String>,
    antisense_genes: Vec<String>,
    diff_to_gene_tss: Option<i32>,
    diff_to_gene_tes: Option<i32>,

    cc: ClassCode,
    ref_gene_id: String,
    ref_gene_name: String,
    ref_tx_id: String,
    ref_length: Option<u32>,
    ref_exon_n: Option<u16>,
    ref_strand: Option<ISOMSTRAND>,
    diff_to_tss: Option<i32>,
    diff_to_tes: Option<i32>,
    query_to_ref_matched_junctions: usize,
    query_to_ref_matched_exons: usize,
    bite: bool,

    all_canonical: bool,
    _rts_stage: bool,
    perc_a_downstream_tts: Option<f32>,
    seq_a_downstream_tts: Option<String>,
    _poly_a_motif: Option<String>,
    _poly_a_motif_found: bool,
    _poly_a_dist: Option<i32>,

    within_cage_peak: Option<bool>,
    dist_to_cage_peak: Option<i32>,
    within_poly_a_peak: Option<bool>,
    dist_to_poly_a_site: Option<i32>,
}

impl ClassifyRecord {
    /// Build a full structural classification for one query transcript.

    pub fn new(
        query_ptir: &QueryPTIR,
        ref_ptir_manager: &RefPTIRManager,
        args: &ClassifyArgs,
    ) -> Self {
        let mut record = Self::empty(query_ptir);

        // Evidence layer 1: reference catalog facts independent of the best transcript hit.
        record.collect_global_reference_evidence(query_ptir, ref_ptir_manager);

        if !record.has_chr_in_ref_gtf {
            return record;
        }

        // chose the primary hit
        let (primary_hit, associated_hits) = find_primary_hit(
            query_ptir,
            ref_ptir_manager,
            &mut record.antisense_genes,
            args,
        );

        if let Some(hit) = primary_hit {
            record.same_strand_overlap_genes = associated_hits
                .iter()
                .map(|hit| hit.gene_id.clone())
                .collect();

            record.apply_primary_hit(&hit, query_ptir);

            record.ref_gene_id = unique_strings(&record.same_strand_overlap_genes).join("_");
            let associated_gene_names: Vec<String> = associated_hits
                .iter()
                .map(|hit| hit.gene_name.clone())
                .collect();
            record.ref_gene_name = unique_strings(&associated_gene_names).join("_");

            record.refine_class_code(query_ptir, ref_ptir_manager);
        } else {
            record.apply_no_same_strand_hit(query_ptir);
        }

        // Final annotation pass: distances to the nearest known gene TSS/TES.
        record.update_gene_endpoint_diffs(query_ptir, ref_ptir_manager);
        record
    }

    fn empty(query_ptir: &QueryPTIR) -> Self {
        Self {
            query_ptir: query_ptir.clone(),
            query_length: query_ptir.transcript_len(),
            query_exon_n: query_ptir.n_exons(),
            query_strand: *query_ptir.strand(),

            has_chr_in_ref_gtf: true,
            per_junction_known_vec: Vec::new(),
            all_junction_known: false,
            per_junction_left_site_known: Vec::new(),
            per_junction_right_site_known: Vec::new(),
            all_splice_sites_known: false,
            contained_in_known_intron: false,
            has_intron_retention_against_catalog: false,
            same_strand_overlap_genes: Vec::new(),
            antisense_genes: Vec::new(),
            diff_to_gene_tss: None,
            diff_to_gene_tes: None,

            cc: ClassCode::Intergenic,
            ref_gene_id: String::new(),
            ref_gene_name: String::new(),
            ref_tx_id: "novel".to_string(),
            ref_length: None,
            ref_exon_n: None,
            ref_strand: None,
            diff_to_tss: None,
            diff_to_tes: None,
            query_to_ref_matched_junctions: 0,
            query_to_ref_matched_exons: 0,
            bite: false,

            all_canonical: matches!(query_ptir.tx_type(), TxType::ALLC | TxType::MONO),
            _rts_stage: false,
            perc_a_downstream_tts: None,
            seq_a_downstream_tts: None,
            _poly_a_motif: None,
            _poly_a_motif_found: false,
            _poly_a_dist: None,

            within_cage_peak: None,
            dist_to_cage_peak: None,
            within_poly_a_peak: None,
            dist_to_poly_a_site: None,
        }
    }

    fn collect_global_reference_evidence(
        &mut self,
        query_ptir: &QueryPTIR,
        ref_ptir_manager: &RefPTIRManager,
    ) {
        if !ref_ptir_manager.has_chr(&query_ptir.chr_name) {
            self.has_chr_in_ref_gtf = false;
            return;
        }

        self.has_chr_in_ref_gtf = true;

        if let Some(junctions) = query_ptir.junction_vec().as_deref() {
            // junction level match
            let (all_junction_known, junction_known_vec) = ref_ptir_manager.junction_match(
                &query_ptir.chr_name,
                junctions,
                query_ptir.strand(),
            );
            self.per_junction_known_vec = junction_known_vec;
            self.all_junction_known = all_junction_known;

            // d/a level match
            let (all_sites_known, left_site_vec, right_site_vec) =
                ref_ptir_manager.splice_site_match(&query_ptir.chr_name, junctions);
            self.per_junction_left_site_known = left_site_vec;
            self.per_junction_right_site_known = right_site_vec;
            self.all_splice_sites_known = all_sites_known;

            // exons cover a junction? -> intron retention against reference transcripts
            self.has_intron_retention_against_catalog = ref_ptir_manager
                .has_intron_retention_against_catalog(
                    &query_ptir.chr_name,
                    query_ptir.strand(),
                    &query_ptir.exons_vec(),
                );
        }

        //  entire transcripts with in a junction -> genic intron
        self.contained_in_known_intron = ref_ptir_manager.contained_in_known_intron(
            &query_ptir.chr_name,
            query_ptir.strand(),
            query_ptir.start(),
            query_ptir.end(),
        );
    }

    /// Copy a selected candidate hit into final output fields.
    fn apply_primary_hit(&mut self, hit: &CandidateHit, _query_ptir: &QueryPTIR) {
        self.cc = match hit.pre_class {
            PreClass::Fsm(sub) => ClassCode::FSM(sub),
            PreClass::Ism(sub) => ClassCode::ISM(sub),
            PreClass::Nic(sub) => ClassCode::NIC(sub),
            PreClass::AnyKnownJunction | PreClass::AnyKnownSpliceSite => {
                ClassCode::NNC(SubNNC::AtLeastOneNovelSpliceSite)
            }
            PreClass::GeneOverlap => ClassCode::Genic,
            PreClass::None => ClassCode::Intergenic,
        };

        self.ref_gene_id = hit.gene_id.clone();
        self.ref_gene_name = hit.gene_name.clone();
        self.ref_tx_id = hit.tx_id.clone();
        self.ref_length = Some(hit.ref_length);
        self.ref_exon_n = Some(hit.ref_exons);
        self.ref_strand = Some(hit.ref_strand);
        self.diff_to_tss = hit.diff_tss;
        self.diff_to_tes = hit.diff_tes;
        self.query_to_ref_matched_junctions = hit.matched_junctions;

        // Exon matches are approximated from matched junctions for splice-based
        // hits; mono/exon-overlap hits report at least one matched exon.
        self.query_to_ref_matched_exons = if self.query_to_ref_matched_junctions > 0 {
            self.query_to_ref_matched_junctions + 1
        } else if hit.exon_overlap_bases > 0 {
            1
        } else {
            0
        };
    }

    /// Refine weak `anyKnown*` evidence into final NIC/NNC/fusion/moreJunctions categories.
    fn refine_class_code(&mut self, query_ptir: &QueryPTIR, ref_ptir_manager: &RefPTIRManager) {
        if !matches!(self.cc, ClassCode::NNC(SubNNC::AtLeastOneNovelSpliceSite)) {
            return;
        }

        let genes = unique_strings(&self.same_strand_overlap_genes);
        if genes.len() == 1 {
            if self.all_splice_sites_known {
                let gene_junctions = ref_ptir_manager
                    .gene_junctions(&genes[0])
                    .unwrap_or_default();
                let gene_junctions: HashSet<(u32, u32)> = gene_junctions.into_iter().collect();
                let all_junctions_in_gene = query_ptir
                    .junction_vec()
                    .as_deref()
                    .unwrap_or(&[])
                    .iter()
                    .all(|junction| gene_junctions.contains(junction));

                self.cc = if self.has_intron_retention_against_catalog {
                    ClassCode::NIC(SubNIC::IntronRetention)
                } else if all_junctions_in_gene {
                    ClassCode::NIC(SubNIC::CombinationOfKnownJunctions)
                } else {
                    ClassCode::NIC(SubNIC::CombinationOfKnownSpliceSites)
                };
            } else {
                self.cc = ClassCode::NNC(SubNNC::AtLeastOneNovelSpliceSite);
            }
        } else if genes.len() > 1 {
            // With multiple associated genes,
            // either be fusion or moreJunctions depending on whether any junction is shared by multiple genes.
            self.cc = if has_junction_shared_by_multiple_genes(query_ptir, &genes, ref_ptir_manager)
            {
                ClassCode::MoreJunctions
            } else {
                ClassCode::Fusion
            };
            self.ref_gene_id = genes.join("_");
            self.ref_gene_name = genes.join("_");
            self.ref_tx_id = "novel".to_string();
            self.ref_length = None;
            self.ref_exon_n = None;
            self.ref_strand = None;
            self.diff_to_tss = None;
            self.diff_to_tes = None;
        }
    }

    /// Resolve transcripts that have no same-strand reference hit.
    /// antisense first, then genic_intron, then intergenic.
    fn apply_no_same_strand_hit(&mut self, query_ptir: &QueryPTIR) {
        if !self.antisense_genes.is_empty() {
            self.cc = ClassCode::Antisense;
            self.same_strand_overlap_genes = unique_strings(&self.antisense_genes)
                .into_iter()
                .map(|gene| format!("novelGene_{}_AS", gene))
                .collect();
            self.ref_gene_id = self.same_strand_overlap_genes.join("_");
            self.ref_gene_name = self.ref_gene_id.clone();
            self.ref_tx_id = "novel".to_string();
        } else if self.contained_in_known_intron {
            self.cc = ClassCode::GenicIntron;
            self.ref_gene_id = "novel".to_string();
            self.ref_gene_name = "novel".to_string();
            self.ref_tx_id = "novel".to_string();
        } else {
            self.cc = ClassCode::Intergenic;
            self.ref_gene_id = "novel".to_string();
            self.ref_gene_name = "novel".to_string();
            self.ref_tx_id = "novel".to_string();
        }

        if query_ptir.n_exons() == 1 && matches!(self.cc, ClassCode::Intergenic) {
            self.ref_length = None;
        }
    }

    /// Fill nearest gene-level TSS/TES distances for associated genes.
    fn update_gene_endpoint_diffs(
        &mut self,
        query_ptir: &QueryPTIR,
        ref_ptir_manager: &RefPTIRManager,
    ) {
        let mut nearest_start: Option<i32> = None;
        let mut nearest_end: Option<i32> = None;

        for gene in &self.same_strand_overlap_genes {
            let Some((starts, ends)) = ref_ptir_manager.gene_starts_ends(gene) else {
                continue;
            };
            for &start in starts {
                update_nearest(&mut nearest_start, start as i32 - query_ptir.start() as i32);
            }
            for &end in ends {
                update_nearest(&mut nearest_end, query_ptir.end() as i32 - end as i32);
            }
        }

        match query_ptir.strand() {
            ISOMSTRAND::Plus | ISOMSTRAND::Unknown => {
                self.diff_to_gene_tss = nearest_start;
                self.diff_to_gene_tes = nearest_end;
            }
            ISOMSTRAND::Minus => {
                self.diff_to_gene_tss = nearest_end;
                self.diff_to_gene_tes = nearest_start;
            }
        }
    }

    /// Write this record as one tab-delimited classification row.
    pub fn write_to_file(&self, writer: &mut dyn Write) -> Result<(), io::Error> {
        let strand = |s: Option<ISOMSTRAND>| {
            s.map(char::from)
                .map(|c| c.to_string())
                .unwrap_or_else(|| "NA".to_string())
        };
        let opt_bool = |b: Option<bool>| match b {
            Some(true) => "true".to_string(),
            Some(false) => "false".to_string(),
            None => "NA".to_string(),
        };
        let opt_f32 = |f: Option<f32>| {
            f.map(|v| format!("{:.4}", v))
                .unwrap_or_else(|| "NA".to_string())
        };
        let opt_i32 = |i: Option<i32>| i.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string());
        let opt_u32 = |i: Option<u32>| i.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string());
        let opt_u16 = |i: Option<u16>| i.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string());
        let opt_str = |s: &Option<String>| s.as_deref().unwrap_or("NA").to_string();

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_ptir.base.source_txid,
            self.query_ptir.chr_name,
            char::from(self.query_strand),
            self.query_length,
            self.query_exon_n,
            self.cc.main_category(),
            self.cc.sub_category(self.query_exon_n),
            empty_as_na(&self.ref_gene_id),
            empty_as_na(&self.ref_gene_name),
            empty_as_na(&self.ref_tx_id),
            opt_u32(self.ref_length),
            opt_u16(self.ref_exon_n),
            strand(self.ref_strand),
            opt_i32(self.diff_to_tss),
            opt_i32(self.diff_to_tes),
            self.query_to_ref_matched_junctions,
            self.query_to_ref_matched_exons,
            self.bite,
            self.all_canonical,
            opt_f32(self.perc_a_downstream_tts),
            opt_str(&self.seq_a_downstream_tts),
            opt_bool(self.within_cage_peak),
            opt_i32(self.dist_to_cage_peak),
            opt_bool(self.within_poly_a_peak),
            opt_i32(self.dist_to_poly_a_site),
        )
    }
}

/// Find the primary same-strand reference hit and all non-overlapping associated hits.

fn find_primary_hit(
    query_ptir: &QueryPTIR,
    ref_ptir_manager: &RefPTIRManager,
    antisense_genes: &mut Vec<String>,
    args: &ClassifyArgs,
) -> (Option<CandidateHit>, Vec<CandidateHit>) {
    let mut best_by_gene: HashMap<String, CandidateHit> = HashMap::default();

    // Gather both mono- and multi-exon reference transcripts whose outer
    // genomic span overlaps the query.
    let refs = ref_ptir_manager.find_overlapping_refs(
        &query_ptir.chr_name,
        query_ptir.start(),
        query_ptir.end(),
    );

    for reference in refs {
        if query_ptir.strand() != reference.standard() {
            // Opposite-strand overlaps are not same-strand hits; they are kept
            // for the later antisense fall-through if no same-strand gene wins.
            antisense_genes.push(reference.base.source_geneid.clone());
            continue;
        }

        // First-stage SQANTI3-like comparison against this one reference transcript.
        let pre_class = classify_against_ref(query_ptir, reference, args);
        if matches!(pre_class, PreClass::None) {
            continue;
        }

        // Convert the pre-class into a scored, metadata-rich candidate hit.
        let hit = CandidateHit::from_ref(query_ptir, reference, pre_class);

        // SQANTI3 keeps the best hit within each gene before comparing genes.
        best_by_gene
            .entry(hit.gene_id.clone())
            .and_modify(|best| {
                if hit.is_better_than(best, query_ptir.n_exons(), query_ptir.transcript_len()) {
                    *best = hit.clone();
                }
            })
            .or_insert(hit);
    }

    let mut hits: Vec<CandidateHit> = best_by_gene.into_values().collect();
    hits.retain(|hit| hit.pre_class.rank() > 0);
    if hits.is_empty() {
        return (None, Vec::new());
    }

    let query_len = query_ptir.transcript_len();
    // Primary ordering is class rank. Composite score is only a tie-breaker.
    hits.sort_by(|a, b| {
        b.pre_class.rank().cmp(&a.pre_class.rank()).then_with(|| {
            b.score(query_ptir.n_exons(), query_len)
                .partial_cmp(&a.score(query_ptir.n_exons(), query_len))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    });

    let mut associated_hits = Vec::new();
    let primary = hits[0].clone();
    let mut current_start = primary.ref_start;
    let mut current_end = primary.ref_end;
    associated_hits.push(primary.clone());

    // Preserve additional gene hits only if their reference spans do not overlap
    // the current associated span, matching SQANTI3's multi-gene association rule.
    for hit in hits.into_iter().skip(1) {
        if hit.pre_class.rank() == 0 {
            break;
        }
        if genomic_overlap_bases(current_start, current_end, hit.ref_start, hit.ref_end) == 0 {
            current_start = current_start.min(hit.ref_start);
            current_end = current_end.max(hit.ref_end);
            associated_hits.push(hit);
        }
    }

    (Some(primary), associated_hits)
}

/// Classify one query against one same-strand reference transcript.
fn classify_against_ref(query: &QueryPTIR, reference: &RefPTIR, args: &ClassifyArgs) -> PreClass {
    let q_mono = query.n_exons() == 1;
    let r_mono = reference.n_exons() == 1;
    let query_exons = query.exons_vec();
    let ref_exons = reference.exons_vec();

    if q_mono && r_mono {
        // Single-exon query against single-exon reference: exonic overlap is the
        // strongest available evidence and is treated as mono-exon FSM.
        return if exon_overlap_bases(&query_exons, &ref_exons) > 0 {
            PreClass::Fsm(SubFSM::MonoExon)
        } else {
            PreClass::None
        };
    }

    if q_mono {
        if exon_overlap_bases(&query_exons, &ref_exons) == 0 {
            return PreClass::None;
        }
        if mono_query_contained_in_ref_exon(query, reference) {
            // Mono-exon query fully contained in a reference exon is an ISM mono-exon case.
            return PreClass::Ism(SubISM::MonoExon);
        }
        if mono_query_spans_ref_intron(query, reference) {
            // Mono-exon query spanning a known intron is NIC intron-retention evidence.
            return PreClass::Nic(SubNIC::MonoExonByIntronRetention);
        }
        return PreClass::GeneOverlap;
    }

    if r_mono {
        // Multi-exon query versus mono-exon reference can only contribute gene overlap.
        return if exon_overlap_bases(&query_exons, &ref_exons) > 0 {
            PreClass::GeneOverlap
        } else {
            PreClass::None
        };
    }

    // Multi-exon query versus multi-exon reference: classify by ordered junction chain.
    match classify_junction_chain(query, reference) {
        JunctionMatch::Exact => {
            let (diff_tss, diff_tes) = same_strand_transcript_space_tss_tes_diffs(query, reference);
            PreClass::Fsm(fsm_subtype(diff_tss, diff_tes, args))
        }
        JunctionMatch::Subset => PreClass::Ism(ism_subtype(query, reference)),
        JunctionMatch::AnyKnownJunction => PreClass::AnyKnownJunction,
        JunctionMatch::AnyKnownSpliceSite => PreClass::AnyKnownSpliceSite,
        JunctionMatch::ExonOverlap => PreClass::GeneOverlap,
        JunctionMatch::NoMatch => PreClass::None,
    }
}

fn fsm_subtype(diff_tss: i32, diff_tes: i32, args: &ClassifyArgs) -> SubFSM {
    // const END_MATCH_BP: i32 = args.fsm_end_match_bp as i32;
    match (
        diff_tss.abs() <= args.fsm_end_match_bp,
        diff_tes.abs() <= args.fsm_end_match_bp,
    ) {
        (true, true) => SubFSM::ReferenceMatch,
        (true, false) => SubFSM::Alternative3End,
        (false, true) => SubFSM::Alternative5End,
        (false, false) => SubFSM::Alternative3And5End,
    }
}

/// decide the ISM subtype
fn ism_subtype(query: &QueryPTIR, reference: &RefPTIR) -> SubISM {
    if query_has_intron_retention_against_ref(query, reference) {
        return SubISM::IntronRetention;
    }

    let q_junctions = query.junction_vec_ref();
    let r_junctions = reference.junction_vec_ref();
    let Some(offset) = find_consecutive_junction_chain(q_junctions, r_junctions) else {
        return SubISM::InternalFragment;
    };

    let agree_front = offset == 0;
    let agree_end = offset + q_junctions.len() == r_junctions.len();
    match (agree_front, agree_end, query.strand()) {
        (true, true, _) => SubISM::Complete,
        (true, false, ISOMSTRAND::Plus | ISOMSTRAND::Unknown) => SubISM::FivePrimeFragment,
        (true, false, ISOMSTRAND::Minus) => SubISM::ThreePrimeFragment,
        (false, true, ISOMSTRAND::Plus | ISOMSTRAND::Unknown) => SubISM::ThreePrimeFragment,
        (false, true, ISOMSTRAND::Minus) => SubISM::FivePrimeFragment,
        (false, false, _) => SubISM::InternalFragment,
    }
}

/// Return true if any query junction is used by more than one associated gene.
/// if a junction from query exactly exists in multiple genes, then its the morejunciton
/// otherwise its a fusion.
fn has_junction_shared_by_multiple_genes(
    query_ptir: &QueryPTIR,
    genes: &[String],
    ref_ptir_manager: &RefPTIRManager,
) -> bool {
    let Some(query_junctions) = query_ptir.junction_vec().as_deref() else {
        return false;
    };

    query_junctions.iter().any(|junction| {
        let mut count = 0;
        for gene in genes {
            let gene_junctions = ref_ptir_manager.gene_junctions(gene).unwrap_or_default();
            if gene_junctions.contains(junction) {
                count += 1;
            }
        }
        count > 1
    })
}

fn update_nearest(current: &mut Option<i32>, value: i32) {
    if current.map_or(true, |existing| value.abs() < existing.abs()) {
        *current = Some(value);
    }
}

fn unique_strings(values: &[String]) -> Vec<String> {
    let mut seen = HashSet::default();
    let mut out = Vec::new();
    for value in values {
        if seen.insert(value.clone()) {
            out.push(value.clone());
        }
    }
    out
}

/// Render empty output identifiers as `NA`.
fn empty_as_na(value: &str) -> &str {
    if value.is_empty() { "NA" } else { value }
}

/// Add sequence-derived QC evidence from the reference FASTA.
pub fn update_group3_seq_context(
    class: &mut ClassifyRecord,
    query_ptir: &QueryPTIR,
    ref_fa: &mut FastaReader,
    args: &ClassifyArgs,
) {
    let chr = &query_ptir.chr_name;
    let strand = query_ptir.strand();
    let tes = query_ptir.base.tes();

    class.all_canonical = matches!(query_ptir.tx_type(), TxType::ALLC | TxType::MONO);

    let chr_len = ref_fa.seq_len(chr).unwrap_or(0);

    let downstream_seq: Option<Vec<u8>> = match strand {
        ISOMSTRAND::Plus => {
            let start = tes as usize;
            let end = start + args.downstream_len;
            (end <= chr_len)
                .then(|| ref_fa.fetch(chr, start, end, false).ok())
                .flatten()
        }
        ISOMSTRAND::Minus => {
            let end = tes as usize;
            (end >= args.downstream_len)
                .then(|| {
                    ref_fa
                        .fetch(chr, end - args.downstream_len, end, false)
                        .ok()
                })
                .flatten()
                .map(|seq| rev_comp(&seq))
        }
        ISOMSTRAND::Unknown => None,
    };

    if let Some(seq) = downstream_seq {
        let a_count = seq
            .iter()
            .filter(|&&base| base == b'A' || base == b'a')
            .count();
        class.perc_a_downstream_tts = Some(a_count as f32 / args.downstream_len as f32 * 100.0);
        class.seq_a_downstream_tts = Some(String::from_utf8_lossy(&seq).into_owned());
    }

    let upstream_seq: Option<Vec<u8>> = match strand {
        ISOMSTRAND::Plus => {
            let end = tes as usize;
            (end >= args.motif_search_window)
                .then(|| {
                    ref_fa
                        .fetch(chr, end - args.motif_search_window, end, false)
                        .ok()
                })
                .flatten()
        }
        ISOMSTRAND::Minus => {
            let start = tes as usize;
            let end = start + args.motif_search_window;
            (end <= chr_len)
                .then(|| ref_fa.fetch(chr, start, end, false).ok())
                .flatten()
                .map(|seq| rev_comp(&seq))
        }
        ISOMSTRAND::Unknown => None,
    };

    if let Some(seq) = upstream_seq {
        let seq_upper: Vec<u8> = seq.iter().map(|base| base.to_ascii_uppercase()).collect();

        for &motif in MOTIFS {
            if let Some(pos) = seq_upper
                .windows(motif.len())
                .position(|window| window == motif)
            {
                class._poly_a_motif = Some(String::from_utf8_lossy(motif).into_owned());
                class._poly_a_motif_found = true;
                class._poly_a_dist =
                    Some((pos + motif.len()) as i32 - args.motif_search_window as i32);
                break;
            }
        }
    }
}

// third parity regions
pub fn update_group4_3rd_party(
    class: &mut ClassifyRecord,
    query_ptir: &QueryPTIR,
    ref_tss: Option<&GuideDb>,
    ref_tes: Option<&GuideDb>,
    _args: &ClassifyArgs,
) {
    let chr = &query_ptir.chr_name;
    let strand = *query_ptir.strand();
    let tss = query_ptir.base.tss();
    let tes = query_ptir.base.tes();

    if let Some(cage_db) = ref_tss {
        // CAGE evidence is evaluated around the query TSS.
        class.within_cage_peak = Some(!cage_db.query_overlaps(chr, strand, tss).is_empty());
        class.dist_to_cage_peak = cage_db
            .query_overlaps_with_flank(chr, &strand, tss, 10_000)
            .into_iter()
            .min_by_key(|iv| {
                let mid = (iv.start as u64 + iv.end as u64) / 2;
                (tss as i64 - mid as i64).unsigned_abs()
            })
            .map(|iv| {
                let mid = ((iv.start as u64 + iv.end as u64) / 2) as i32;
                let dist = tss as i32 - mid;
                if strand == ISOMSTRAND::Minus {
                    -dist
                } else {
                    dist
                }
            });
    }

    if let Some(polya_db) = ref_tes {
        // PolyA peak evidence is evaluated around the query TES/TTS.
        class.within_poly_a_peak = Some(!polya_db.query_overlaps(chr, strand, tes).is_empty());
        class.dist_to_poly_a_site = polya_db
            .query_overlaps_with_flank(chr, &strand, tes, 10_000)
            .into_iter()
            .min_by_key(|iv| {
                let mid = (iv.start as u64 + iv.end as u64) / 2;
                (tes as i64 - mid as i64).unsigned_abs()
            })
            .map(|iv| {
                let mid = ((iv.start as u64 + iv.end as u64) / 2) as i32;
                let dist = tes as i32 - mid;
                if strand == ISOMSTRAND::Minus {
                    -dist
                } else {
                    dist
                }
            });
    }
}
