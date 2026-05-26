// use std::collections::HashSet;

use std::io::{self, Write};

use ahash::HashSet;

use crate::{
    classify::{class_code::*, query_ptir::QueryPTIR, ref_ptir::RefPTIR},
    core::{tx_strand::ISOMSTRAND, tx_type::TxType},
    index::fasta::FastaReader,
    merge::guide::GuideDb,
    utils::rev_comp,
};
#[derive(Debug, Clone)]
pub struct ClassifyRecord {
    // Group 0: query PTIR
    query_ptir: QueryPTIR,

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
    pub fn new(query_ptir: &QueryPTIR) -> Self {
        Self {
            query_ptir: query_ptir.clone(),

            query_length: query_ptir.end() - query_ptir.start(),
            query_exon_n: query_ptir.n_exons() as u32,
            query_strand: *query_ptir.standard(),

            cc: ClassCode::Intergenic,
            ref_gene_id: String::new(),
            ref_gene_name: String::new(),
            ref_tx_id: String::new(),
            ref_length: 0,
            ref_exon_n: 0,
            ref_strand: ISOMSTRAND::Unknown,
            diff_to_tss: 0,
            diff_to_tes: 0,
            _diff_to_gene_tss: 0,
            _diff_to_gene_tts: 0,
            query_to_ref_matched_junctions: 0,
            query_to_ref_matched_exons: 0,
            bite: false,

            all_canonical: false,
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

    pub fn new_intergenic(query_ptir: &QueryPTIR) -> Self {
        Self {
            query_ptir: query_ptir.clone(),

            query_length: query_ptir.end() - query_ptir.start(),
            query_exon_n: query_ptir.n_exons() as u32,
            query_strand: *query_ptir.standard(),

            cc: ClassCode::Intergenic,
            ref_gene_id: String::new(),
            ref_gene_name: String::new(),
            ref_tx_id: String::new(),
            ref_length: 0,
            ref_exon_n: 0,
            ref_strand: ISOMSTRAND::Unknown,
            diff_to_tss: 0,
            diff_to_tes: 0,
            _diff_to_gene_tss: 0,
            _diff_to_gene_tts: 0,
            query_to_ref_matched_junctions: 0,
            query_to_ref_matched_exons: 0,
            bite: false,

            all_canonical: false,
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

    pub fn new_fusion(query_ptir: &QueryPTIR, gene_set: &HashSet<&str>) -> Self {
        let gene_names: Vec<String> = gene_set.iter().map(|x| x.to_string()).collect();

        let gene_names_string = gene_names.join("_");

        Self {
            query_ptir: query_ptir.clone(),

            query_length: query_ptir.end() - query_ptir.start(),
            query_exon_n: query_ptir.n_exons() as u32,
            query_strand: *query_ptir.standard(),

            cc: ClassCode::Fusion,
            ref_gene_id: "NA".to_string(),
            ref_gene_name: gene_names_string,
            ref_tx_id: "novel".to_string(),
            ref_length: 0,
            ref_exon_n: 0,
            ref_strand: ISOMSTRAND::Unknown,
            diff_to_tss: 0,
            diff_to_tes: 0,
            _diff_to_gene_tss: 0,
            _diff_to_gene_tts: 0,
            query_to_ref_matched_junctions: 0,
            query_to_ref_matched_exons: 0,
            bite: false,

            all_canonical: false,
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

    pub fn write_to_file(&self, writer: &mut dyn Write) -> Result<(), io::Error> {
        let sc = |s: ISOMSTRAND| char::from(s);
        let opt_bool = |b: Option<bool>| {
            if b == Some(true) {
                "true"
            } else if b == Some(false) {
                "false"
            } else {
                "NA"
            }
        };
        let opt_f32 = |f: Option<f32>| {
            f.map(|v| format!("{:.4}", v))
                .unwrap_or_else(|| "NA".to_string())
        };
        let opt_str = |s: &Option<String>| s.as_deref().unwrap_or("NA").to_string();
        let opt_i32 = |i: Option<i32>| i.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string());
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_ptir.base.source_txid,
            self.query_ptir.chr_name,
            sc(self.query_strand),
            self.query_length,
            self.query_exon_n,
            self.cc,
            self.cc.sub_category(),
            self.ref_gene_id,
            self.ref_gene_name,
            self.ref_tx_id,
            self.ref_length,
            self.ref_exon_n,
            sc(self.ref_strand),
            self.diff_to_tss,
            self.diff_to_tes,
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

impl ClassifyRecord {
    pub(super) fn cc(&self) -> ClassCode {
        self.cc
    }

    pub(super) fn endpoint_total_diff(&self) -> i32 {
        self.diff_to_tss.abs() + self.diff_to_tes.abs()
    }

    pub(super) fn matched_junctions(&self) -> usize {
        self.query_to_ref_matched_junctions
    }

    pub(super) fn exon_count_diff(&self) -> u32 {
        (self.query_exon_n as i32 - self.ref_exon_n as i32).abs() as u32
    }

    pub fn ref_gene_name(&self) -> &str {
        &self.ref_gene_name
    }
}

fn find_consecutive_subseq(query: &[(u32, u32)], reference: &[(u32, u32)]) -> Option<usize> {
    reference.windows(query.len()).position(|w| w == query)
}

pub fn classify(
    query_ptir: &QueryPTIR,
    ref_ptir: &RefPTIR,
    // ref_tss: Option<&GuideDb>,
    // ref_tes: Option<&GuideDb>,
    // ref_fa: &mut FastaReader,
) -> ClassifyRecord {
    let mut class = ClassifyRecord::new(query_ptir);

    class.cc = get_class_code(query_ptir, ref_ptir);

    // Group 2
    class.ref_tx_id = ref_ptir.base.source_txid.clone();
    class.ref_gene_id = ref_ptir.base.source_geneid.clone();
    class.ref_gene_name = ref_ptir.gene_name.clone();
    class.ref_length = ref_ptir.end() - ref_ptir.start();
    class.ref_exon_n = ref_ptir.n_exons() as u32;
    class.ref_strand = *ref_ptir.standard();
    class.diff_to_tss = query_ptir.base.tss() as i32 - ref_ptir.base.tss() as i32; // plus number, larger, on right side.
    class.diff_to_tes = query_ptir.base.tes() as i32 - ref_ptir.base.tes() as i32;
    //
    class.bite = query_ptir.start() > ref_ptir.start() && query_ptir.end() < ref_ptir.end();

    if let (Some(qj), Some(rj)) = (query_ptir.junction_vec(), ref_ptir.junction_vec()) {
        let ref_set: HashSet<(u32, u32)> = rj.iter().copied().collect();
        class.query_to_ref_matched_junctions = qj.iter().filter(|j| ref_set.contains(j)).count();
        if class.query_to_ref_matched_junctions > 0 {
            class.query_to_ref_matched_exons = class.query_to_ref_matched_junctions + 1;
        }
    } else {
        // at least one is mono-exon; overlapping counts as 1
        class.query_to_ref_matched_exons = 1;
    }

    class
}

pub fn get_class_code(query_ptir: &QueryPTIR, ref_ptir: &RefPTIR) -> ClassCode {
    let q_strand = query_ptir.standard();
    let r_strand = ref_ptir.standard();

    // Antisense: opposite known strands
    if *q_strand != ISOMSTRAND::Unknown && *r_strand != ISOMSTRAND::Unknown && q_strand != r_strand
    {
        return ClassCode::Antisense;
    }

    let q_mono = query_ptir.n_exons() == 1;
    let r_mono = ref_ptir.n_exons() == 1;

    // Both mono-exon
    if q_mono && r_mono {
        return ClassCode::FSM(SubFSM::MonoExon);
    }

    // Mono query vs multi-exon ref
    if q_mono {
        let ref_junctions = ref_ptir.junction_vec().as_deref().unwrap();
        let q_start = query_ptir.start();
        let q_end = query_ptir.end();

        // NIC: query spans a complete ref intron
        for (intron_s, intron_e) in ref_junctions {
            if q_start <= *intron_s && q_end >= *intron_e {
                return ClassCode::NIC(SubNIC::MonoExonIntronRetention);
            }
        }

        // ISM: query is fully within a single ref exon
        let mut exon_start = ref_ptir.start();
        for (intron_s, intron_e) in ref_junctions {
            if q_start >= exon_start && q_end <= *intron_s {
                return ClassCode::ISM(SubISM::MonoExon);
            }
            exon_start = *intron_e;
        }
        if q_start >= exon_start && q_end <= ref_ptir.end() {
            return ClassCode::ISM(SubISM::MonoExon);
        }

        // Determine whether query touches exonic and/or intronic regions.
        // NIC already ruled out complete intron spanning, so intron overlaps here are partial.
        let mut touches_exon = false;
        let mut touches_intron = false;
        let mut exon_start = ref_ptir.start();
        for (intron_s, intron_e) in ref_junctions {
            if q_start < *intron_s && q_end > exon_start {
                touches_exon = true;
            }
            if q_start < *intron_e && q_end > *intron_s {
                touches_intron = true;
            }
            exon_start = *intron_e;
        }
        if q_start < ref_ptir.end() && q_end > exon_start {
            touches_exon = true;
        }

        return match (touches_exon, touches_intron) {
            (false, true) => ClassCode::GenicIntron,
            (true, true) => ClassCode::GenicGenomic,
            _ => ClassCode::NNC(SubNNC::MonoExon),
        };
    }

    // Multi-exon query vs mono-exon ref → no known splice sites
    if r_mono {
        return ClassCode::NNC(SubNNC::AtLeastOneNovelSpliceSite);
    }

    // Both multi-exon
    let q_junctions = query_ptir.junction_vec().as_deref().unwrap();
    let r_junctions = ref_ptir.junction_vec().as_deref().unwrap();

    // FSM: identical junction chains
    if q_junctions == r_junctions {
        const THRESHOLD: i32 = 50;
        let diff_tss = (query_ptir.base.tss() as i32 - ref_ptir.base.tss() as i32).abs();
        let diff_tes = (query_ptir.base.tes() as i32 - ref_ptir.base.tes() as i32).abs();
        return match (diff_tss <= THRESHOLD, diff_tes <= THRESHOLD) {
            (true, true) => ClassCode::FSM(SubFSM::ReferenceMatch),
            (true, false) => ClassCode::FSM(SubFSM::Alternative3End),
            (false, true) => ClassCode::FSM(SubFSM::Alternative5End),
            (false, false) => ClassCode::FSM(SubFSM::Alternative5And3End),
        };
    }

    // ISM: query junctions are a strict consecutive subsequence of ref
    if q_junctions.len() < r_junctions.len() {
        if let Some(offset) = find_consecutive_subseq(q_junctions, r_junctions) {
            let is_5_complete = offset == 0;
            let is_3_complete = offset + q_junctions.len() == r_junctions.len();

            // Check if a ref intron is retained at the fragment boundary
            let mut retention = false;
            if !is_3_complete {
                let next_intron = r_junctions[offset + q_junctions.len()];
                if query_ptir.end() > next_intron.0 {
                    retention = true;
                }
            }
            if !retention && !is_5_complete {
                let prev_intron = r_junctions[offset - 1];
                if query_ptir.start() < prev_intron.1 {
                    retention = true;
                }
            }
            if retention {
                return ClassCode::ISM(SubISM::IntronRetention);
            }
            return match (is_5_complete, is_3_complete) {
                (true, false) => ClassCode::ISM(SubISM::FivePrimeFragment),
                (false, true) => ClassCode::ISM(SubISM::ThreePrimeFragment),
                _ => ClassCode::ISM(SubISM::InternalFragment),
            };
        }
    }

    // NIC / NNC: build ref site sets
    let ref_donors: HashSet<u32> = r_junctions.iter().map(|(s, _)| *s).collect();
    let ref_acceptors: HashSet<u32> = r_junctions.iter().map(|(_, e)| *e).collect();
    let ref_junc_set: HashSet<(u32, u32)> = r_junctions.iter().copied().collect();

    let mut all_donors_known = true;
    let mut all_acceptors_known = true;
    let mut all_pairs_known = true;
    for (d, a) in q_junctions {
        if !ref_donors.contains(d) {
            all_donors_known = false;
        }
        if !ref_acceptors.contains(a) {
            all_acceptors_known = false;
        }
        if !ref_junc_set.contains(&(*d, *a)) {
            all_pairs_known = false;
        }
    }

    if all_donors_known && all_acceptors_known {
        // Check for intron retention: a query exon fully spans a ref intron
        let mut q_exon_start = query_ptir.start();
        let has_retention = {
            let mut found = false;
            for (qi_s, qi_e) in q_junctions {
                let q_exon_end = *qi_s;
                if r_junctions
                    .iter()
                    .any(|(ri_s, ri_e)| *ri_s >= q_exon_start && *ri_e <= q_exon_end)
                {
                    found = true;
                    break;
                }
                q_exon_start = *qi_e;
            }
            found
        };
        if has_retention {
            return ClassCode::NIC(SubNIC::IntronRetention);
        }
        if all_pairs_known {
            return ClassCode::NIC(SubNIC::CombinationOfKnownJunctions);
        }
        return ClassCode::NIC(SubNIC::CombinationOfKnownSpliceSites);
    }

    ClassCode::NNC(SubNNC::AtLeastOneNovelSpliceSite)
}

pub fn update_group3(class: &mut ClassifyRecord, query_ptir: &QueryPTIR, ref_fa: &mut FastaReader) {
    let chr = &query_ptir.chr_name;
    let strand = query_ptir.standard();
    let tts = query_ptir.base.tes();

    class.all_canonical = matches!(query_ptir.tx_type(), TxType::ALLC);

    let chr_len = ref_fa.seq_len(chr).unwrap_or(0);

    // perc_a_downstream_tts + seq_a_downstream_tts: 20bp downstream of TTS, strand-corrected
    const DOWNSTREAM_LEN: usize = 20;
    let downstream_seq: Option<Vec<u8>> = match strand {
        // make this assert in fetch
        ISOMSTRAND::Plus => {
            let start = tts as usize;
            let end = start + DOWNSTREAM_LEN;
            if end <= chr_len {
                ref_fa.fetch(chr, start, end, false).ok()
            } else {
                None
            }
        }
        ISOMSTRAND::Minus => {
            let end = tts as usize;
            if end >= DOWNSTREAM_LEN {
                ref_fa
                    .fetch(chr, end - DOWNSTREAM_LEN, end, false)
                    .ok()
                    .map(|s| rev_comp(&s))
            } else {
                None
            }
        }
        ISOMSTRAND::Unknown => None,
    };
    if let Some(seq) = downstream_seq {
        let a_count = seq.iter().filter(|&&b| b == b'A' || b == b'a').count();
        class.perc_a_downstream_tts = Some(a_count as f32 / DOWNSTREAM_LEN as f32 * 100.0);
        class.seq_a_downstream_tts = Some(String::from_utf8_lossy(&seq).into_owned());
    }

    // _rts_stage: check for direct repeat ≥4bp at each junction donor boundary
    if let Some(junctions) = query_ptir.junction_vec() {
        const RTS_WINDOW: usize = 25;
        const MIN_REPEAT: usize = 4;
        class._rts_stage = junctions.iter().any(|(intron_start, _)| {
            let donor = *intron_start as usize;
            if donor < RTS_WINDOW || donor + RTS_WINDOW > chr_len {
                return false;
            }
            let window = ref_fa
                .fetch(chr, donor - RTS_WINDOW, donor + RTS_WINDOW, false)
                .unwrap_or_default();
            let (exon_tail, intron_head) = window.split_at(RTS_WINDOW.min(window.len()));
            (MIN_REPEAT..=RTS_WINDOW).any(|k| {
                k <= exon_tail.len()
                    && k <= intron_head.len()
                    && exon_tail[exon_tail.len() - k..].eq_ignore_ascii_case(&intron_head[..k])
            })
        });
    }

    // poly_a motif: search the 50bp upstream of TTS (strand-corrected)
    const MOTIF_WINDOW: usize = 50;
    let upstream_seq: Option<Vec<u8>> = match strand {
        ISOMSTRAND::Plus => {
            let end = tts as usize;
            if end >= MOTIF_WINDOW {
                ref_fa.fetch(chr, end - MOTIF_WINDOW, end, false).ok()
            } else {
                None
            }
        }
        ISOMSTRAND::Minus => {
            let start = tts as usize;
            let end = start + MOTIF_WINDOW;
            if end <= chr_len {
                ref_fa
                    .fetch(chr, start, end, false)
                    .ok()
                    .map(|s| rev_comp(&s))
            } else {
                None
            }
        }
        ISOMSTRAND::Unknown => None,
    };
    if let Some(seq) = upstream_seq {
        let seq_upper: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
        const MOTIFS: &[&[u8]] = &[
            b"AATAAA", b"ATTAAA", b"AGTAAA", b"TATAAA", b"CATAAA", b"GATAAA", b"AATATA", b"AATACA",
            b"AATAGA", b"AAAAAG", b"ACTAAA", b"AAGAAA", b"AATGAA", b"TTTAAA", b"AAAACA", b"GGGGCT",
        ];
        for &motif in MOTIFS {
            if let Some(pos) = seq_upper.windows(motif.len()).position(|w| w == motif) {
                let dist = (pos + motif.len()) as i32 - MOTIF_WINDOW as i32;
                class._poly_a_motif = Some(String::from_utf8_lossy(motif).into_owned());
                class._poly_a_motif_found = true;
                class._poly_a_dist = Some(dist);
                break;
            }
        }
    }
}

pub fn update_group4(
    class: &mut ClassifyRecord,
    query_ptir: &QueryPTIR,
    ref_tss: Option<&GuideDb>,
    ref_tes: Option<&GuideDb>,
) {
    let chr = &query_ptir.chr_name;
    let strand = *query_ptir.standard();
    let tss = query_ptir.base.tss();
    let tes = query_ptir.base.tes();

    if let Some(cage_db) = ref_tss {
        class.within_cage_peak = Some(!cage_db.query_overlaps(chr, strand, tss).is_empty());
        let near_peaks = cage_db.query_overlaps_with_flank(chr, &strand, tss, 10_000);
        class.dist_to_cage_peak = near_peaks
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
        class.within_poly_a_peak = Some(!polya_db.query_overlaps(chr, strand, tes).is_empty());
        let near_sites = polya_db.query_overlaps_with_flank(chr, &strand, tes, 10_000);
        class.dist_to_poly_a_site = near_sites
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
