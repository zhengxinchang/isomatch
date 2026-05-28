use core::fmt;

use crate::{
    MergeArgs,
    core::{ptir::PTIR, tx_strand::ISOMSTRAND, tx_type::TxType},
    merge::{
        grouped_ptirs::{GroupedPTIR, GroupedPTIREntry},
        guide::GuideDb,
        merge_error::MergeError,
    },
};
use clap::ValueEnum;
use rustc_hash::FxHashMap;
use serde::Serialize;

#[derive(Copy, Clone, Debug, Serialize, ValueEnum)]
pub enum MergePolicyArg {
    Longer,
    Shorter,
    Major,
}

#[derive(Copy, Clone, Debug, Serialize)]
pub enum MergePolicyUsed {
    Longer,
    Shorter,
    Major,
    Guide(GuideResolution),
}

impl fmt::Display for MergePolicyUsed {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            MergePolicyUsed::Longer => "longer",
            MergePolicyUsed::Shorter => "shorter",
            MergePolicyUsed::Major => "major",
            MergePolicyUsed::Guide(resolution) => match resolution {
                GuideResolution::Definitive => "guide_definitive",
                GuideResolution::Majority => "guide_majority",
                GuideResolution::Longer => "guide_longer",
            },
        };
        write!(f, "{s}")
    }
}

impl MergePolicyUsed {
    pub fn from_arg_policy(arg_policy: &MergePolicyArg) -> Self {
        match *arg_policy {
            MergePolicyArg::Longer => MergePolicyUsed::Longer,
            MergePolicyArg::Shorter => MergePolicyUsed::Shorter,
            MergePolicyArg::Major => MergePolicyUsed::Major,
        }
    }
}

#[derive(Copy, Clone, Debug, Serialize)]
pub enum GuideResolution {
    Definitive, // 所有 guide 候选位置一致
    Majority,   // guide 过滤后多数投票决定
    Longer,     // guide 过滤后仍平票，用 Longer 兜底
}

#[derive(Copy, Clone, Debug, Serialize, ValueEnum)]
pub enum TerminalRefineMode {
    None,
    TSS,
    TES,
    Both,
}

pub fn merge_cluster(
    chrom: &str,
    n_exon: u16,
    strand: ISOMSTRAND,
    cluster_idx: &Vec<usize>,
    scluster: &[PTIR],
    args: &MergeArgs,
    guide_tss: &Option<GuideDb>,
    guide_tes: &Option<GuideDb>,
) -> Result<Vec<GroupedPTIR>, MergeError> {
    // match strand {
    //     ISOMSTRAND::Plus | ISOMSTRAND::Minus | ISOMSTRAND::Unknown => {
    if n_exon != 1 {
        // multiple exon tx
        // further splice ptirs into canonical and non-canonical
        let mut canonical_global_idxs = Vec::new();
        let mut non_canonical_global_idxs = Vec::new();
        for ptir_idx in cluster_idx {
            match scluster[*ptir_idx].tx_type {
                TxType::ALLC => {
                    canonical_global_idxs.push(ptir_idx);
                }
                TxType::NOTC | TxType::PRTC => {
                    non_canonical_global_idxs.push(*ptir_idx);
                }
                TxType::MONO => {
                    return Err(MergeError::TxType {
                        reason: "Mono exon transcript shouldn't have exon number > 1".to_string(),
                    });
                }
            }
        }
        // merge canonical

        // println!("canonical_vec:{:?}", &canonical_global_idxs,);

        let grpptirs: Vec<GroupedPTIR> =
            merge_canonical(canonical_global_idxs, scluster, args, &strand, n_exon)?;

        // furhter split grpptirs based on tss and tes in args
        let mut grpptirs: Vec<GroupedPTIR> = refine_canonical_grouped_ptir(grpptirs, &strand, args);

        // process the canonical
        for grpptir in grpptirs.iter_mut() {
            grpptir.profile_canonical_ptirs(chrom, args, guide_tss, guide_tes)?;
        }

        // println!("noncanonical_vec:{:?}", &non_canonical_global_idxs,);
        let rest_non_canonical_global_idxs = noncannonical_to_canonical(
            &mut grpptirs,
            scluster,
            non_canonical_global_idxs,
            &strand,
            args,
        )?;

        let rest_grpptirs = match rest_non_canonical_global_idxs {
            Some(ptir_idxs) => {
                // println!("rest_grpptirs:{:?}", &ptir_idxs);
                merge_rest_noncanonical(ptir_idxs, scluster, &strand, args)?
            }
            None => Vec::new(),
        };

        let mut rest_grpptirs = refine_non_canonical_grouped_ptir(rest_grpptirs, &strand, args);

        for grpptir in rest_grpptirs.iter_mut() {
            grpptir.profile_non_canonical_ptirs(chrom, args, guide_tss, guide_tes)?;
        }

        grpptirs.extend(rest_grpptirs.into_iter());
        return Ok(grpptirs);
    } else {
        return merge_mono_exon(
            chrom,
            cluster_idx,
            scluster,
            &strand,
            args,
            guide_tss,
            guide_tes,
        );
    }
}

// Merge canonical transcripts into at least one MPTIR, each one is a
// backbone
// for canonical transcirpt, no wobble are allowed by default.
// union find algorithm to find allow wobble match of SJ.
pub fn merge_canonical(
    canonical_vec: Vec<&usize>,
    super_cluster: &[PTIR],
    args: &MergeArgs,
    strand: &ISOMSTRAND,
    n_exon: u16,
) -> Result<Vec<GroupedPTIR>, MergeError> {
    // the ptir has same strand, same exons and all canonical
    // only consider the TSS and TES distance is ok
    // because all canonical should be able to correctly splice

    // in case of the no canonical txs
    if canonical_vec.is_empty() {
        return Ok(Vec::new());
    }

    let mut sorted_tx_indice: Vec<usize> = canonical_vec.iter().map(|&&i| i).collect();

    // sort tx indice by the left most junction.
    sorted_tx_indice.sort_by(|&a, &b| {
        let ja = super_cluster[a].junction_vec.as_ref().unwrap();
        let jb = super_cluster[b].junction_vec.as_ref().unwrap();
        if ja[0].0 == jb[0].0 {
            if ja[0].1 == jb[0].1 {
                std::cmp::Ordering::Equal
            } else {
                ja[0].1.cmp(&jb[0].1)
            }
        } else {
            ja[0].0.cmp(&jb[0].0)
        }
    });

    let mut grpptirs = Vec::new();

    let n = sorted_tx_indice.len();
    let mut parent: Vec<usize> = (0..n).collect();
    let mut size = vec![1usize; n];

    for left in 0..n {
        let left_junc = super_cluster[sorted_tx_indice[left]]
            .junctions()
            .ok_or(MergeError::NoJunctionFound)?;

        for right in left + 1..n {
            let right_junc = super_cluster[sorted_tx_indice[right]]
                .junctions()
                .ok_or(MergeError::NoJunctionFound)?;
            let (in_wobble, _, _) = is_splice_junctions_match(
                left_junc, right_junc, strand, args.wob_a, args.wob_d, args.wob_u,
            );
            if in_wobble {
                uf_union(&mut parent, &mut size, left, right);
            }
        }
    }

    let mut groups: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for local_idx in 0..n {
        let root = uf_find(&mut parent, local_idx);
        groups.entry(root).or_default().push(local_idx);
    }

    let mut grouped_local_indices: Vec<Vec<usize>> = groups.into_values().collect();
    grouped_local_indices.sort_by_key(|group| group[0]);

    for group in grouped_local_indices {
        let mut grpptir = GroupedPTIR::new(strand, n_exon);

        for &local_idx in group.iter() {
            let super_idx = sorted_tx_indice[local_idx];
            grpptir.add_canonical_ptir(&super_cluster[super_idx], super_idx)?;
        }

        grpptirs.push(grpptir);
    }

    Ok(grpptirs)
}

pub fn refine_canonical_grouped_ptir(
    grouped_ptirs: Vec<GroupedPTIR>,
    _strand: &ISOMSTRAND,
    args: &MergeArgs,
) -> Vec<GroupedPTIR> {
    if matches!(args.terminal_refine, TerminalRefineMode::None) {
        return grouped_ptirs;
    }

    let mut refined = Vec::new();
    for grpptir in grouped_ptirs {
        let strand = grpptir.strand();
        let n_exon = grpptir.n_exon();
        let entries = grpptir.canonical_entries_cloned();
        let split_entries = split_grouped_entries_by_terminals(
            entries,
            strand,
            args.tss_wob,
            args.tes_wob,
            args.terminal_refine,
        );

        for entries in split_entries {
            refined.push(GroupedPTIR::from_canonical_entries(
                &strand, n_exon, entries,
            ));
        }
    }
    refined
}

pub fn refine_non_canonical_grouped_ptir(
    grouped_ptirs: Vec<GroupedPTIR>,
    _strand: &ISOMSTRAND,
    args: &MergeArgs,
) -> Vec<GroupedPTIR> {
    if matches!(args.terminal_refine_nc, TerminalRefineMode::None) {
        return grouped_ptirs;
    }

    let mut refined = Vec::new();
    for grpptir in grouped_ptirs {
        let strand = grpptir.strand();
        let n_exon = grpptir.n_exon();
        let entries = grpptir.non_canonical_entries_cloned();
        let split_entries = split_grouped_entries_by_terminals(
            entries,
            strand,
            args.tss_wob_nc,
            args.tes_wob_nc,
            args.terminal_refine_nc,
        );

        for entries in split_entries {
            refined.push(GroupedPTIR::from_non_canonical_entries(
                &strand, n_exon, entries,
            ));
        }
    }
    refined
}

fn split_grouped_entries_by_terminals(
    mut entries: Vec<GroupedPTIREntry>,
    strand: ISOMSTRAND,
    tss_wob: u32,
    tes_wob: u32,
    mode: TerminalRefineMode,
) -> Vec<Vec<GroupedPTIREntry>> {
    if entries.len() <= 1 || matches!(mode, TerminalRefineMode::None) {
        return vec![entries];
    }

    entries.sort_by_key(|entry| {
        let (tss, tes) = entry_terminals(entry.left, entry.right, strand);
        (tss, tes)
    });

    let mut groups: Vec<Vec<GroupedPTIREntry>> = Vec::new();
    let mut current_group: Vec<GroupedPTIREntry> = Vec::new();
    let mut anchor_tss = 0;
    let mut anchor_tes = 0;

    for entry in entries {
        let (curr_tss, curr_tes) = entry_terminals(entry.left, entry.right, strand);

        if current_group.is_empty() {
            anchor_tss = curr_tss;
            anchor_tes = curr_tes;
            current_group.push(entry);
            continue;
        }

        if terminal_match_with_anchor(
            curr_tss, curr_tes, anchor_tss, anchor_tes, tss_wob, tes_wob, mode,
        ) {
            current_group.push(entry);
        } else {
            groups.push(current_group);
            current_group = vec![entry];
            anchor_tss = curr_tss;
            anchor_tes = curr_tes;
        }
    }

    if !current_group.is_empty() {
        groups.push(current_group);
    }

    groups
}

fn entry_terminals(left: u32, right: u32, strand: ISOMSTRAND) -> (u32, u32) {
    match strand {
        ISOMSTRAND::Plus => (left, right),
        ISOMSTRAND::Minus => (right, left),
        ISOMSTRAND::Unknown => (left, right),
    }
}

fn terminal_match_with_anchor(
    curr_tss: u32,
    curr_tes: u32,
    anchor_tss: u32,
    anchor_tes: u32,
    tss_wob: u32,
    tes_wob: u32,
    mode: TerminalRefineMode,
) -> bool {
    match mode {
        TerminalRefineMode::None => true,
        TerminalRefineMode::TSS => curr_tss.abs_diff(anchor_tss) <= tss_wob,
        TerminalRefineMode::TES => curr_tes.abs_diff(anchor_tes) <= tes_wob,
        TerminalRefineMode::Both => {
            curr_tss.abs_diff(anchor_tss) <= tss_wob && curr_tes.abs_diff(anchor_tes) <= tes_wob
        }
    }
}

pub fn noncannonical_to_canonical(
    grpptirs: &mut [GroupedPTIR],
    super_cluster: &[PTIR],
    noncanonical_vec: Vec<usize>,
    strand: &ISOMSTRAND,
    args: &MergeArgs,
) -> Result<Option<Vec<usize>>, MergeError> {
    if grpptirs.is_empty() {
        return Ok(Some(noncanonical_vec));
    }

    // only select the best grpptir for each no-all-ca ptir
    let mut marked_ptirs = FxHashMap::default();
    noncanonical_vec.iter().for_each(|idx| {
        marked_ptirs.insert(*idx, (false, usize::MAX, u32::MAX, u32::MAX)); // absorbed, grpptirs_local_idx, d_diff_bp,a_diff_bp
    });

    for &nidx in noncanonical_vec.iter() {
        let ptir = &super_cluster[nidx];
        // let ptir_juncs = get_junctions(ptir)?;
        let ptir_juncs = ptir.junctions().ok_or(MergeError::NoJunctionFound)?;
        for (grpptirs_local_idx, grpptir) in grpptirs.iter_mut().enumerate() {
            let mptir_junction = grpptir.repr_junction();
            match is_splice_junctions_match(
                ptir_juncs,
                mptir_junction,
                strand,
                args.wob_a_nc,
                args.wob_d_nc,
                args.wob_u_nc,
            ) {
                (true, new_d_diff_bp, new_a_diff_bp) => {
                    // grpptir.add_non_canonical_ptir(ptir, nidx)?;

                    marked_ptirs.entry(nidx).and_modify(
                        |(flag, grp_local_idx, d_diff_bp, a_diff_bp)| {
                            if should_replace_noncanonical_match(
                                *flag,
                                new_d_diff_bp,
                                new_a_diff_bp,
                                *d_diff_bp,
                                *a_diff_bp,
                            ) {
                                *flag = true;
                                *grp_local_idx = grpptirs_local_idx;
                                *d_diff_bp = new_d_diff_bp;
                                *a_diff_bp = new_a_diff_bp;
                            }
                        },
                    );
                }
                (false, _, _) => {
                    // rest_unmerged_ptirs.push(nidx);
                }
            }
        }
    }
    let mut rest = Vec::new();
    for (ptir_idx, (flag, grptir_idx, _, _)) in marked_ptirs.into_iter() {
        if flag == true {
            // absorbton
            grpptirs[grptir_idx].add_non_canonical_ptir(&super_cluster[ptir_idx], ptir_idx)?;
        } else {
            rest.push(ptir_idx);
        }
    }

    rest.dedup_by(|a, b| a == b);
    Ok(Some(rest))
}

fn should_replace_noncanonical_match(
    has_current_match: bool,
    new_d_diff_bp: u32,
    new_a_diff_bp: u32,
    current_d_diff_bp: u32,
    current_a_diff_bp: u32,
) -> bool {
    if !has_current_match {
        return true;
    }

    let new_total = u64::from(new_d_diff_bp) + u64::from(new_a_diff_bp);
    let current_total = u64::from(current_d_diff_bp) + u64::from(current_a_diff_bp);
    new_total < current_total
}

pub fn merge_rest_noncanonical(
    rest_non_canonical_ptirs: Vec<usize>,
    scluster: &[PTIR],
    strand: &ISOMSTRAND,
    args: &MergeArgs,
) -> Result<Vec<GroupedPTIR>, MergeError> {
    if rest_non_canonical_ptirs.is_empty() {
        return Ok(Vec::new());
    }

    let mut sorted_tx_indice: Vec<usize> = rest_non_canonical_ptirs.iter().copied().collect();

    sorted_tx_indice.sort_by(|&a, &b| {
        let ja = scluster[a].junction_vec.as_ref().unwrap();
        let jb = scluster[b].junction_vec.as_ref().unwrap();
        if ja[0].0 == jb[0].0 {
            if ja[0].1 == jb[0].1 {
                std::cmp::Ordering::Equal
            } else {
                ja[0].1.cmp(&jb[0].1)
            }
        } else {
            ja[0].0.cmp(&jb[0].0)
        }
    });

    let n = sorted_tx_indice.len();
    let mut parent: Vec<usize> = (0..n).collect();
    let mut size = vec![1usize; n];

    for left in 0..n {
        let left_junc = scluster[sorted_tx_indice[left]]
            .junctions()
            .ok_or(MergeError::NoJunctionFound)?;

        for right in left + 1..n {
            let right_junc = scluster[sorted_tx_indice[right]]
                .junctions()
                .ok_or(MergeError::NoJunctionFound)?;
            let (in_wobble, _, _) = is_splice_junctions_match(
                left_junc,
                right_junc,
                strand,
                args.wob_a_nc,
                args.wob_d_nc,
                args.wob_u_nc,
            );
            if in_wobble {
                uf_union(&mut parent, &mut size, left, right);
            }
        }
    }

    let mut groups: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for local_idx in 0..n {
        let root = uf_find(&mut parent, local_idx);
        groups.entry(root).or_default().push(local_idx);
    }

    let mut grouped_local_indices: Vec<Vec<usize>> = groups.into_values().collect();
    grouped_local_indices.sort_by_key(|group| group[0]);

    let first_super_idx = sorted_tx_indice[0];
    let n_exon = scluster[first_super_idx].n_exons;
    let mut grpptirs = Vec::new();

    for group in grouped_local_indices {
        let mut grpptir = GroupedPTIR::new(strand, n_exon);

        for &local_idx in &group {
            let super_idx = sorted_tx_indice[local_idx];
            grpptir.add_non_canonical_ptir(&scluster[super_idx], super_idx)?;
        }

        grpptirs.push(grpptir);
    }

    Ok(grpptirs)
}

pub fn merge_mono_exon(
    chrom: &str,
    scluster_idxs: &[usize],
    scluster: &[PTIR],
    strand: &ISOMSTRAND,
    args: &MergeArgs,
    guide_tss: &Option<GuideDb>,
    guide_tes: &Option<GuideDb>,
) -> Result<Vec<GroupedPTIR>, MergeError> {
    // 对scluster_idxs 按照PTIR的start 和end 排序
    if scluster_idxs.is_empty() {
        return Ok(Vec::new());
    }

    let mut sorted_tx_indice: Vec<usize> = scluster_idxs.to_vec();
    sorted_tx_indice.sort_by_key(|&idx| (scluster[idx].start, scluster[idx].end));

    let n = sorted_tx_indice.len();
    let mut parent: Vec<usize> = (0..n).collect();
    let mut size = vec![1usize; n];

    for left in 0..n {
        let left_ptir = &scluster[sorted_tx_indice[left]];
        for right in left + 1..n {
            let right_ptir = &scluster[sorted_tx_indice[right]];

            if right_ptir.start > left_ptir.end {
                break;
            }

            if mono_reciprocal_overlap(left_ptir, right_ptir) >= args.mono_ovlp {
                uf_union(&mut parent, &mut size, left, right);
            }
        }
    }

    let mut groups: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for local_idx in 0..n {
        let root = uf_find(&mut parent, local_idx);
        groups.entry(root).or_default().push(local_idx);
    }

    let mut grouped_local_indices: Vec<Vec<usize>> = groups.into_values().collect();
    grouped_local_indices.sort_by_key(|group| group[0]);

    let mut grpptirs = Vec::new();
    for group in grouped_local_indices {
        let mut grpptir = GroupedPTIR::new(strand, 1);
        for &local_idx in &group {
            let super_idx = sorted_tx_indice[local_idx];
            grpptir.add_mono_ptir(&scluster[super_idx], super_idx);
        }
        grpptir.profile_mono_ptirs(chrom, args, guide_tss, guide_tes)?;
        grpptirs.push(grpptir);
    }

    Ok(grpptirs)
}

fn mono_reciprocal_overlap(left: &PTIR, right: &PTIR) -> f64 {
    let overlap_start = left.start.max(right.start);
    let overlap_end = left.end.min(right.end);

    if overlap_start > overlap_end {
        return 0.0;
    }

    let overlap_len = (overlap_end - overlap_start + 1) as f64;
    let left_len = (left.end - left.start + 1) as f64;
    let right_len = (right.end - right.start + 1) as f64;

    let left_frac = overlap_len / left_len;
    let right_frac = overlap_len / right_len;

    left_frac.min(right_frac)
}

fn uf_find(parent: &mut [usize], node: usize) -> usize {
    if parent[node] != node {
        let root = uf_find(parent, parent[node]);
        parent[node] = root;
    }
    parent[node]
}

fn uf_union(parent: &mut [usize], size: &mut [usize], left: usize, right: usize) {
    let mut left_root = uf_find(parent, left);
    let mut right_root = uf_find(parent, right);

    if left_root == right_root {
        return;
    }

    if size[left_root] < size[right_root] {
        std::mem::swap(&mut left_root, &mut right_root);
    }

    parent[right_root] = left_root;
    size[left_root] += size[right_root];
}

pub fn is_splice_junctions_match(
    curr: &[(u32, u32)],
    other: &[(u32, u32)],
    strand: &ISOMSTRAND,
    // args: &MergeArgs,
    awob: u32,
    dwob: u32,
    uwob: u32,
) -> (
    bool, // if in wobble under current parameters
    u32,  // how many bp difference in donor
    u32,  //how many bp difference in acceptor
) {
    if curr.len() != other.len() {
        return (false, u32::MAX, u32::MAX);
    }

    let mut bp_diff_a = 0;
    let mut bp_diff_d = 0;
    let mut in_wobble = true;
    curr.iter()
        .zip(other.iter())
        .for_each(|(curr_junc, other_junc)| match strand {
            ISOMSTRAND::Plus => {
                let donor_diff = curr_junc.0.abs_diff(other_junc.0);
                let acceptor_diff = curr_junc.1.abs_diff(other_junc.1);
                bp_diff_a += acceptor_diff;
                bp_diff_d += donor_diff;
                if donor_diff > dwob || acceptor_diff > awob {
                    in_wobble = false
                }
            }
            ISOMSTRAND::Minus => {
                let donor_diff = curr_junc.1.abs_diff(other_junc.1);
                let acceptor_diff = curr_junc.0.abs_diff(other_junc.0);
                bp_diff_a += acceptor_diff;
                bp_diff_d += donor_diff;
                if donor_diff > dwob || acceptor_diff > awob {
                    in_wobble = false
                }
            }
            ISOMSTRAND::Unknown => {
                let left_diff = curr_junc.0.abs_diff(other_junc.0);
                let right_diff = curr_junc.1.abs_diff(other_junc.1);
                // For unknown strand we keep a stable left/right convention:
                // left coordinate contributes to donor_diff, right to acceptor_diff.
                bp_diff_d += left_diff;
                bp_diff_a += right_diff;
                if left_diff > uwob || right_diff > uwob {
                    in_wobble = false
                }
            }
        });
    (in_wobble, bp_diff_d, bp_diff_a)
}

#[cfg(test)]
mod tests {
    use super::{is_splice_junctions_match, should_replace_noncanonical_match};
    use crate::core::tx_strand::ISOMSTRAND;

    #[test]
    fn first_noncanonical_match_is_accepted_without_overflow() {
        assert!(should_replace_noncanonical_match(
            false,
            1,
            2,
            u32::MAX,
            u32::MAX,
        ));
    }

    #[test]
    fn worse_noncanonical_match_is_not_selected() {
        assert!(!should_replace_noncanonical_match(true, 5, 5, 1, 1));
    }

    #[test]
    fn unknown_strand_wobble_uses_same_donor_acceptor_convention_as_output_stats() {
        let curr = vec![(100, 200)];
        let other = vec![(103, 207)];
        let (matched, donor_diff, acceptor_diff) =
            is_splice_junctions_match(&curr, &other, &ISOMSTRAND::Unknown, 10, 10, 10);

        assert!(matched);
        assert_eq!(donor_diff, 3);
        assert_eq!(acceptor_diff, 7);
    }
}
