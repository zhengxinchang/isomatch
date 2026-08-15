use super::grouped_ptirs::GroupedPTIR;
use crate::MergeArgs;
use crate::core::tx_strand::ISOMSTRAND;
use rustc_hash::{FxHashMap, FxHashSet};
pub fn assign_global_ids(
    all_grouped_ptirs: &mut [GroupedPTIR],
    global_gene_id: &mut u32,
    global_tx_id: &mut u32,
    args: &MergeArgs,
) {
    all_grouped_ptirs.iter_mut().for_each(|grp| {
        *global_tx_id += 1;
        grp.update_tx_id(*global_tx_id);
    });

    all_grouped_ptirs.sort_by_key(|grp| (grp.strand(), grp.start(), grp.end()));

    for strand in [ISOMSTRAND::Plus, ISOMSTRAND::Minus, ISOMSTRAND::Unknown] {
        let strand_idxs = strand_indices(all_grouped_ptirs, strand);
        for overlap_group in split_by_overlap(all_grouped_ptirs, &strand_idxs) {
            let (connected_groups, unassigned_idxs) =
                split_by_splice_site_connectivity(all_grouped_ptirs, &overlap_group);
            let mut gene_groups = Vec::new();
            let mut bridge_groups = Vec::new();

            for connected_group in connected_groups {
                let (split_gene_groups, bridge_group) =
                    split_read_through_group(all_grouped_ptirs, &connected_group);
                gene_groups.extend(split_gene_groups);
                if !bridge_group.is_empty() {
                    bridge_groups.push(bridge_group);
                }
            }

            let mut remaining_unassigned = assign_unassigned_to_gene_groups(
                all_grouped_ptirs,
                &mut gene_groups,
                &unassigned_idxs,
                args,
            );
            remaining_unassigned
                .sort_by_key(|&idx| (all_grouped_ptirs[idx].start(), all_grouped_ptirs[idx].end()));

            gene_groups.extend(group_unassigned_remains(
                all_grouped_ptirs,
                &remaining_unassigned,
                args,
            ));
            gene_groups.extend(bridge_groups);

            for gene_group in gene_groups {
                update_next_gene_id(all_grouped_ptirs, &gene_group, global_gene_id);
            }
        }
    }
}

fn strand_indices(grps: &[GroupedPTIR], strand: ISOMSTRAND) -> Vec<usize> {
    grps.iter()
        .enumerate()
        .filter_map(|(idx, grp)| (grp.strand() == strand).then_some(idx))
        .collect()
}

fn split_by_overlap(grps: &[GroupedPTIR], idxs: &[usize]) -> Vec<Vec<usize>> {
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut current_group = Vec::new();
    let mut current_max_end = 0;

    for &idx in idxs {
        let grp = &grps[idx];
        if current_group.is_empty() || grp.start() <= current_max_end {
            current_group.push(idx);
            current_max_end = current_max_end.max(grp.end());
        } else {
            groups.push(std::mem::take(&mut current_group));
            current_group.push(idx);
            current_max_end = grp.end();
        }
    }

    if !current_group.is_empty() {
        groups.push(current_group);
    }

    groups
}

type ConnectedGroups = Vec<Vec<usize>>;
type UnassignedIdxs = Vec<usize>;
type Exon = (u32, u32);

fn split_by_splice_site_connectivity(
    grps: &[GroupedPTIR],
    idxs: &[usize],
) -> (ConnectedGroups, UnassignedIdxs) {
    let splice_sites: Vec<FxHashSet<u32>> = idxs
        .iter()
        .map(|&idx| {
            grps[idx]
                .repr_junction()
                .iter()
                .flat_map(|&(left, right)| [left, right])
                .collect()
        })
        .collect();

    let mut parent: Vec<usize> = (0..idxs.len()).collect();
    let mut size = vec![1usize; idxs.len()];
    let mut has_splice_match = vec![false; idxs.len()];

    // each transcirpt's splice sites
    for left in 0..idxs.len() {
        // skip if mono exon
        if splice_sites[left].is_empty() {
            continue;
        }

        // skip if mono exon for rest of the transcirpt
        for right in left + 1..idxs.len() {
            if splice_sites[right].is_empty() {
                continue;
            }

            if splice_sites[left]
                .iter()
                .any(|site| splice_sites[right].contains(site))
            // if any overlap of two transcirpts
            {
                has_splice_match[left] = true;
                has_splice_match[right] = true;
                idx_union(&mut parent, &mut size, left, right);
            }
        }
    }

    let mut root_to_group_idx: FxHashMap<usize, usize> = FxHashMap::default();
    let mut connected_groups: Vec<Vec<usize>> = Vec::new();
    let mut unassigned_idxs = Vec::new();
    for local_idx in 0..idxs.len() {
        if !has_splice_match[local_idx] {
            unassigned_idxs.push(idxs[local_idx]);
            continue;
        }

        let root = idx_find(&mut parent, local_idx);
        if let Some(&group_idx) = root_to_group_idx.get(&root) {
            connected_groups[group_idx].push(idxs[local_idx]);
        } else {
            root_to_group_idx.insert(root, connected_groups.len());
            connected_groups.push(vec![idxs[local_idx]]);
        }
    }

    (connected_groups, unassigned_idxs)
}

fn assign_unassigned_to_gene_groups(
    grps: &[GroupedPTIR],
    gene_groups: &mut [Vec<usize>],
    unassigned_idxs: &[usize],
    args: &MergeArgs,
) -> Vec<usize> {
    let mut gene_group_exons = Vec::with_capacity(gene_groups.len());
    for gene_group in gene_groups.iter() {
        let exons: Vec<Exon> = gene_group
            .iter()
            .flat_map(|&idx| grps[idx].exons_from_repr())
            .collect();
        gene_group_exons.push(merge_exons(exons));
    }

    let mut remaining_unassigned = Vec::new();
    for &unassigned_idx in unassigned_idxs {
        let transcript_exons = grps[unassigned_idx].exons_from_repr();
        let transcript_length = exon_len(&transcript_exons);

        let mut best_group_idx = None;
        let mut best_overlap = 0u64;
        let mut best_is_tied = false;

        for (group_idx, group_exons) in gene_group_exons.iter().enumerate() {
            let overlap = exon_overlap_len(&transcript_exons, group_exons);

            if overlap * 100 < transcript_length * args.gene_assign_min_overlap {
                continue;
            }

            if overlap > best_overlap {
                best_group_idx = Some(group_idx);
                best_overlap = overlap;
                best_is_tied = false;
            } else if overlap == best_overlap {
                best_is_tied = true;
            }
        }

        if let Some(group_idx) = best_group_idx
            && !best_is_tied
        {
            gene_groups[group_idx].push(unassigned_idx);
        } else {
            remaining_unassigned.push(unassigned_idx);
        }
    }

    remaining_unassigned
}

fn group_unassigned_remains(
    grps: &[GroupedPTIR],
    idxs: &[usize],
    args: &MergeArgs,
) -> Vec<Vec<usize>> {
    let transcript_exons: Vec<Vec<(u32, u32)>> = idxs
        .iter()
        .map(|&idx| grps[idx].exons_from_repr())
        .collect();
    let transcript_lengths: Vec<u64> = transcript_exons
        .iter()
        .map(|exons| exon_len(exons))
        .collect();

    let mut parent: Vec<usize> = (0..idxs.len()).collect();
    let mut size = vec![1usize; idxs.len()];

    // The reciprocal-overlap relation is transitive through union-find and can chain groups.
    for left in 0..idxs.len() {
        for right in left + 1..idxs.len() {
            let overlap = exon_overlap_len(&transcript_exons[left], &transcript_exons[right]);

            let passes_left =
                overlap * 100 >= transcript_lengths[left] * args.unassigned_group_min_overlap;
            let passes_right =
                overlap * 100 >= transcript_lengths[right] * args.unassigned_group_min_overlap;
            if passes_left && passes_right {
                idx_union(&mut parent, &mut size, left, right);
            }
        }
    }

    let mut root_to_group_idx: FxHashMap<usize, usize> = FxHashMap::default();
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for (local_idx, &idx) in idxs.iter().enumerate() {
        let root = idx_find(&mut parent, local_idx);
        if let Some(&group_idx) = root_to_group_idx.get(&root) {
            groups[group_idx].push(idx);
        } else {
            root_to_group_idx.insert(root, groups.len());
            groups.push(vec![idx]);
        }
    }

    groups
}

fn merge_exons(mut exons: Vec<Exon>) -> Vec<Exon> {
    exons.sort_unstable_by_key(|&(start, end)| (start, end));

    let mut merged_exons: Vec<Exon> = Vec::new();
    for (start, end) in exons {
        if let Some(last) = merged_exons.last_mut()
            && start <= last.1.saturating_add(1)
        {
            last.1 = last.1.max(end);
        } else {
            merged_exons.push((start, end));
        }
    }

    merged_exons
}

fn exon_len(exons: &[Exon]) -> u64 {
    exons
        .iter()
        .map(|&(start, end)| u64::from(end) - u64::from(start) + 1)
        .sum()
}

fn exon_overlap_len(left: &[Exon], right: &[Exon]) -> u64 {
    let mut overlap = 0u64;
    let mut left_exon_idx = 0;
    let mut right_exon_idx = 0;

    while left_exon_idx < left.len() && right_exon_idx < right.len() {
        let (left_start, left_end) = left[left_exon_idx];
        let (right_start, right_end) = right[right_exon_idx];

        let overlap_start = left_start.max(right_start);
        let overlap_end = left_end.min(right_end);
        if overlap_start <= overlap_end {
            overlap += u64::from(overlap_end) - u64::from(overlap_start) + 1;
        }

        if left_end <= right_end {
            left_exon_idx += 1;
        }
        if right_end <= left_end {
            right_exon_idx += 1;
        }
    }

    overlap
}

fn update_next_gene_id(grps: &mut [GroupedPTIR], idxs: &[usize], curr_gene_id: &mut u32) {
    *curr_gene_id += 1;
    let gene_id = *curr_gene_id;

    for &idx in idxs {
        grps[idx].update_gene_id(gene_id);
    }
}

fn idx_find(parent: &mut [usize], node: usize) -> usize {
    if parent[node] != node {
        let root = idx_find(parent, parent[node]);
        parent[node] = root;
    }
    parent[node]
}

fn idx_union(parent: &mut [usize], size: &mut [usize], left: usize, right: usize) {
    let mut left_root = idx_find(parent, left);
    let mut right_root = idx_find(parent, right);

    if left_root == right_root {
        return;
    }

    // always let left root larger
    if size[left_root] < size[right_root] {
        std::mem::swap(&mut left_root, &mut right_root);
    }

    parent[right_root] = left_root;
    size[left_root] += size[right_root];
}

struct CoverageSegment {
    pos: u64,
    next_pos: u64,
    segment_cov: usize,
}

pub fn split_read_through_group(
    grps: &[GroupedPTIR],
    idxs: &[usize],
) -> (Vec<Vec<usize>>, Vec<usize>) {
    if idxs.len() < 3 {
        return (vec![idxs.to_vec()], Vec::new());
    }

    let profile = build_span_coverage_profile(grps, idxs);
    if profile.len() < 3 {
        return (vec![idxs.to_vec()], Vec::new());
    }

    let Some(valley_idx) = best_read_through_valley(&profile) else {
        return (vec![idxs.to_vec()], Vec::new());
    };

    let cut =
        profile[valley_idx].pos + (profile[valley_idx].next_pos - profile[valley_idx].pos) / 2;
    let (left_group, right_group, bridge_group) = partition_by_cut(grps, idxs, cut);

    if left_group.is_empty() || right_group.is_empty() {
        return (vec![idxs.to_vec()], Vec::new());
    }

    (vec![left_group, right_group], bridge_group)
}

fn build_span_coverage_profile(grps: &[GroupedPTIR], idxs: &[usize]) -> Vec<CoverageSegment> {
    let mut events = Vec::with_capacity(idxs.len() * 2);
    for &idx in idxs {
        let grp = &grps[idx];
        events.push((u64::from(grp.start()), 1i32));
        events.push((u64::from(grp.end()) + 1, -1i32));
    }
    events.sort_unstable_by_key(|event| event.0);

    // Each profile item covers [start, end) with a constant transcript-span coverage.
    let mut profile: Vec<CoverageSegment> = Vec::new();
    let mut coverage = 0i32;
    let mut event_idx = 0;
    while event_idx < events.len() {
        let position = events[event_idx].0;
        while event_idx < events.len() && events[event_idx].0 == position {
            coverage += events[event_idx].1;
            event_idx += 1;
        }

        if event_idx == events.len() {
            break;
        }

        let next_position = events[event_idx].0;
        let segment_coverage = coverage as usize;
        if let Some(last) = profile.last_mut()
            && last.next_pos == position
            && last.segment_cov == segment_coverage
        {
            last.next_pos = next_position;
        } else {
            profile.push(CoverageSegment {
                pos: position,
                next_pos: next_position,
                segment_cov: segment_coverage,
            });
        }
    }

    profile
}

fn best_read_through_valley(profile: &[CoverageSegment]) -> Option<usize> {
    let mut right_peaks = vec![0usize; profile.len()];
    for profile_idx in (0..profile.len() - 1).rev() {
        right_peaks[profile_idx] =
            right_peaks[profile_idx + 1].max(profile[profile_idx + 1].segment_cov);
    }

    let mut left_peak = profile[0].segment_cov;
    let mut best_valley_idx = None;
    let mut best_valley_coverage = usize::MAX;
    let mut best_valley_width = 0;
    for profile_idx in 1..profile.len() - 1 {
        let valley_coverage = profile[profile_idx].segment_cov;
        let right_peak = right_peaks[profile_idx];
        let flank_coverage = left_peak.min(right_peak);

        let is_local_minimum = valley_coverage < profile[profile_idx - 1].segment_cov
            && valley_coverage < profile[profile_idx + 1].segment_cov;
        let is_deep_enough =
            flank_coverage >= 2 && valley_coverage.saturating_mul(2) <= flank_coverage;

        if is_local_minimum && is_deep_enough {
            let valley_width = profile[profile_idx].next_pos - profile[profile_idx].pos;
            if valley_coverage < best_valley_coverage
                || (valley_coverage == best_valley_coverage && valley_width > best_valley_width)
            {
                best_valley_idx = Some(profile_idx);
                best_valley_coverage = valley_coverage;
                best_valley_width = valley_width;
            }
        }

        left_peak = left_peak.max(valley_coverage);
    }

    best_valley_idx
}

fn partition_by_cut(
    grps: &[GroupedPTIR],
    idxs: &[usize],
    cut: u64,
) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let mut left_group = Vec::new();
    let mut right_group = Vec::new();
    let mut bridge_group = Vec::new();
    for &idx in idxs {
        let grp = &grps[idx];
        if u64::from(grp.end()) < cut {
            left_group.push(idx);
        } else if u64::from(grp.start()) >= cut {
            right_group.push(idx);
        } else {
            bridge_group.push(idx);
        }
    }

    (left_group, right_group, bridge_group)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exon_helpers_merge_and_count_inclusive_overlap() {
        let exons = merge_exons(vec![(20, 25), (1, 10), (11, 15)]);
        assert_eq!(exons, vec![(1, 15), (20, 25)]);
        assert_eq!(exon_len(&exons), 21);
        assert_eq!(exon_overlap_len(&exons, &[(5, 12), (18, 22)]), 11);
    }

    #[test]
    fn best_read_through_valley_prefers_deeper_then_wider() {
        let profile = vec![
            CoverageSegment {
                pos: 1,
                next_pos: 10,
                segment_cov: 4,
            },
            CoverageSegment {
                pos: 10,
                next_pos: 20,
                segment_cov: 1,
            },
            CoverageSegment {
                pos: 20,
                next_pos: 30,
                segment_cov: 4,
            },
            CoverageSegment {
                pos: 30,
                next_pos: 50,
                segment_cov: 1,
            },
            CoverageSegment {
                pos: 50,
                next_pos: 55,
                segment_cov: 4,
            },
        ];

        assert_eq!(best_read_through_valley(&profile), Some(3));
    }
}
