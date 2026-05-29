use std::collections::HashSet;

use crate::{
    classify::{query_ptir::QueryPTIR, ref_ptir::RefPTIR},
    core::tx_strand::ISOMSTRAND,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JunctionMatch {
    Exact,
    Subset,
    AnyKnownJunction,
    AnyKnownSpliceSite,
    ExonOverlap,
    NoMatch,
}

pub fn exon_overlap_bases(query: &[(u32, u32)], reference: &[(u32, u32)]) -> u32 {
    let mut total = 0;
    for &(qs, qe) in query {
        for &(rs, re) in reference {
            let start = qs.max(rs);
            let end = qe.min(re);
            if start < end {
                total += end - start;
            }
        }
    }
    total
}

pub fn genomic_overlap_bases(q_start: u32, q_end: u32, r_start: u32, r_end: u32) -> u32 {
    let start = q_start.max(r_start);
    let end = q_end.min(r_end);
    end.saturating_sub(start)
}

pub fn splice_site_agreement(query_exons: &[(u32, u32)], ref_exons: &[(u32, u32)]) -> usize {
    let q_sites: HashSet<u32> = query_exons
        .iter()
        .flat_map(|&(start, end)| [start, end])
        .collect();

    ref_exons
        .iter()
        .flat_map(|&(start, end)| [start, end])
        .filter(|site| q_sites.contains(site))
        .count()
}

pub fn find_consecutive_junction_chain<T: PartialEq>(
    query: &[T],
    reference: &[T],
) -> Option<usize> {
    if query.is_empty() || query.len() > reference.len() {
        return None;
    }
    reference
        .windows(query.len())
        .position(|window| window == query)
}

pub fn classify_junction_chain(query: &QueryPTIR, reference: &RefPTIR) -> JunctionMatch {
    let query_exons = query.exons_vec();
    let ref_exons = reference.exons_vec();

    let Some(q_junctions) = query.junction_vec().as_deref() else {
        return if exon_overlap_bases(&query_exons, &ref_exons) > 0 {
            JunctionMatch::ExonOverlap
        } else {
            JunctionMatch::NoMatch
        };
    };
    let Some(r_junctions) = reference.junction_vec().as_deref() else {
        return if exon_overlap_bases(&query_exons, &ref_exons) > 0 {
            JunctionMatch::ExonOverlap
        } else {
            JunctionMatch::NoMatch
        };
    };

    if q_junctions == r_junctions {
        return JunctionMatch::Exact;
    }

    if q_junctions.len() < r_junctions.len()
        && find_consecutive_junction_chain(q_junctions, r_junctions).is_some()
    {
        return JunctionMatch::Subset;
    }

    let ref_junctions: HashSet<(u32, u32)> = r_junctions.iter().copied().collect();
    if q_junctions
        .iter()
        .any(|junction| ref_junctions.contains(junction))
    {
        return JunctionMatch::AnyKnownJunction;
    }

    if splice_site_agreement(&query_exons, &ref_exons) > 0 {
        return JunctionMatch::AnyKnownSpliceSite;
    }

    if exon_overlap_bases(&query_exons, &ref_exons) > 0 {
        return JunctionMatch::ExonOverlap;
    }

    JunctionMatch::NoMatch
}

// calculate the difference of tss and tes at trasncript space, not genomic space.
// meaning the intron length wont be counted.
pub fn same_strand_transcript_space_tss_tes_diffs(
    query: &QueryPTIR,
    reference: &RefPTIR,
) -> (i32, i32) {
    let min_start = query.start().min(reference.start());
    let relative_q_start = query.start() - min_start;
    let relative_q_end = relative_q_start + query.transcript_len();
    let relative_r_start = reference.start() - min_start;
    let relative_r_end = relative_r_start + reference.transcript_len();

    match query.strand() {
        ISOMSTRAND::Plus | ISOMSTRAND::Unknown => (
            relative_r_start as i32 - relative_q_start as i32,
            relative_q_end as i32 - relative_r_end as i32,
        ),
        ISOMSTRAND::Minus => (
            relative_q_end as i32 - relative_r_end as i32,
            relative_r_start as i32 - relative_q_start as i32,
        ),
    }
}

pub fn mono_query_contained_in_ref_exon(query: &QueryPTIR, reference: &RefPTIR) -> bool {
    let q_start = query.start();
    let q_end = query.end();
    reference
        .exons_vec()
        .iter()
        .any(|&(start, end)| start <= q_start && q_end <= end)
}

pub fn mono_query_spans_ref_intron(query: &QueryPTIR, reference: &RefPTIR) -> bool {
    let q_start = query.start();
    let q_end = query.end();
    reference
        .junction_vec_ref()
        .iter()
        .any(|&(start, end)| q_start < start && start < end && end < q_end)
}

pub fn query_has_intron_retention_against_ref(query: &QueryPTIR, reference: &RefPTIR) -> bool {
    let q_exons = query.exons_vec();
    let r_junctions = reference.junction_vec_ref();
    q_exons.iter().any(|&(q_start, q_end)| {
        r_junctions
            .iter()
            .any(|&(r_start, r_end)| q_start <= r_start && r_start < r_end && r_end <= q_end)
    })
}
