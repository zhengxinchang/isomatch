use std::io::Write;

use crate::{
    MergeArgs,
    core::{ptir::PTIR, status::TxType, tx_strand::ISOMSTRAND},
    merge::{merge_error::MergeError, policy::MergePolicy},
};
use rustc_hash::FxHashMap;

pub(crate) type GroupedPTIREntry = (
    usize,           // super cluster idx
    u32,             // left （start）
    u32,             // right （end）
    Vec<(u32, u32)>, // junction
);

pub struct GroupedPTIR {
    // status:MPTIRTYPE,
    gene_id: u32,
    tx_id: u32,
    strand: ISOMSTRAND,
    n_exon: u16,
    all_canonical_ptir_counts: u32,
    all_canonical_ptir_list: Vec<GroupedPTIREntry>,
    canonical_junction_range: Vec<(u32, u32, u32, u32)>, // the range of each junction among all canonical ptirs, (min_left, max_left, min_right, max_right)
    no_all_canonical_ptir_counts: u32,
    no_all_canonical_ptir_list: Vec<GroupedPTIREntry>,
    repr_junction: Vec<(u32, u32)>,
    repr_left: u32,
    repr_right: u32,
    used_repr_junction_policy: MergePolicy,
    used_repr_left_policy: MergePolicy,
    used_repr_right_policy: MergePolicy,
    used_repr_mono_policy: MergePolicy,
    repr_loaded: bool,
}

impl GroupedPTIR {
    pub fn new(strand: &ISOMSTRAND, n_exon: u16) -> GroupedPTIR {
        GroupedPTIR {
            gene_id: 0,
            tx_id: 0,
            strand: strand.clone(),
            n_exon: n_exon,
            all_canonical_ptir_counts: 0,
            all_canonical_ptir_list: Vec::new(),
            canonical_junction_range: Vec::new(),
            no_all_canonical_ptir_counts: 0,
            no_all_canonical_ptir_list: Vec::new(),
            repr_junction: Vec::new(),
            repr_left: 0,
            repr_right: 0,
            used_repr_junction_policy: MergePolicy::Major,
            used_repr_left_policy: MergePolicy::Major,
            used_repr_right_policy: MergePolicy::Major,
            used_repr_mono_policy: MergePolicy::Major,
            repr_loaded: false,
        }
    }

    pub fn repr_tss(&self) -> u32 {
        boundaries_to_terminals(self.repr_left, self.repr_right, self.strand).0
    }

    pub fn repr_tes(&self) -> u32 {
        boundaries_to_terminals(self.repr_left, self.repr_right, self.strand).1
    }

    pub fn tss(&self) -> u32 {
        self.repr_tss()
    }

    pub fn tes(&self) -> u32 {
        self.repr_tes()
    }

    pub(crate) fn strand(&self) -> ISOMSTRAND {
        self.strand
    }

    pub(crate) fn n_exon(&self) -> u16 {
        self.n_exon
    }

    pub(crate) fn canonical_entries_cloned(&self) -> Vec<GroupedPTIREntry> {
        self.all_canonical_ptir_list.clone()
    }

    pub(crate) fn non_canonical_entries_cloned(&self) -> Vec<GroupedPTIREntry> {
        self.no_all_canonical_ptir_list.clone()
    }

    pub(crate) fn from_canonical_entries(
        strand: &ISOMSTRAND,
        n_exon: u16,
        entries: Vec<GroupedPTIREntry>,
    ) -> GroupedPTIR {
        GroupedPTIR {
            gene_id: 0,
            tx_id: 0,
            strand: *strand,
            n_exon,
            all_canonical_ptir_counts: entries.len() as u32,
            all_canonical_ptir_list: entries,
            canonical_junction_range: Vec::new(),
            no_all_canonical_ptir_counts: 0,
            no_all_canonical_ptir_list: Vec::new(),
            repr_junction: Vec::new(),
            repr_left: 0,
            repr_right: 0,
            used_repr_junction_policy: MergePolicy::Major,
            used_repr_left_policy: MergePolicy::Major,
            used_repr_right_policy: MergePolicy::Major,
            used_repr_mono_policy: MergePolicy::Major,
            repr_loaded: false,
        }
    }

    pub(crate) fn from_non_canonical_entries(
        strand: &ISOMSTRAND,
        n_exon: u16,
        entries: Vec<GroupedPTIREntry>,
    ) -> GroupedPTIR {
        GroupedPTIR {
            gene_id: 0,
            tx_id: 0,
            strand: *strand,
            n_exon,
            all_canonical_ptir_counts: 0,
            all_canonical_ptir_list: Vec::new(),
            canonical_junction_range: Vec::new(),
            no_all_canonical_ptir_counts: entries.len() as u32,
            no_all_canonical_ptir_list: entries,
            repr_junction: Vec::new(),
            repr_left: 0,
            repr_right: 0,
            used_repr_junction_policy: MergePolicy::Major,
            used_repr_left_policy: MergePolicy::Major,
            used_repr_right_policy: MergePolicy::Major,
            used_repr_mono_policy: MergePolicy::Major,
            repr_loaded: false,
        }
    }

    fn set_repr_from_terminals(&mut self, tss: u32, tes: u32) {
        let (repr_left, repr_right) = terminals_to_boundaries(tss, tes, self.strand);
        self.repr_left = repr_left;
        self.repr_right = repr_right;
    }

    fn set_used_repr_terminal_policies(
        &mut self,
        tss_policy: MergePolicy,
        tes_policy: MergePolicy,
    ) {
        if tss_is_left_boundary(self.strand) {
            self.used_repr_left_policy = tss_policy;
            self.used_repr_right_policy = tes_policy;
        } else {
            self.used_repr_left_policy = tes_policy;
            self.used_repr_right_policy = tss_policy;
        }
    }

    pub fn add_canonical_ptir(
        &mut self,
        ptir: &PTIR,
        scluster_idx: usize,
    ) -> Result<(), MergeError> {
        self.all_canonical_ptir_counts += 1;
        let junction = ptir
            .junction_vec
            .clone()
            .ok_or(MergeError::NoJunctionFound)?;
        self.all_canonical_ptir_list
            .push((scluster_idx, ptir.start, ptir.end, junction));
        Ok(())
    }

    pub fn profile_canonical_ptirs(&mut self, args: &MergeArgs) -> Result<(), MergeError> {
        // generate canonical_junction_range based on canonical_ptir_list
        self.canonical_junction_range.clear();
        self.repr_junction.clear();
        self.repr_left = 0;
        self.repr_right = 0;
        self.used_repr_junction_policy = args.splice_policy;
        self.used_repr_left_policy = args.tss_policy;
        self.used_repr_right_policy = args.tes_policy;
        self.used_repr_mono_policy = args.mono_policy;

        let Some((_, _, _, first_junctions)) = self.all_canonical_ptir_list.first() else {
            return Err(MergeError::NoJunctionFound);
        };

        let junction_count = first_junctions.len();
        self.canonical_junction_range = first_junctions
            .iter()
            .map(|&(left, right)| (left, left, right, right))
            .collect();

        for (_, _, _, junctions) in self.all_canonical_ptir_list.iter().skip(1) {
            debug_assert_eq!(
                junctions.len(),
                self.canonical_junction_range.len(),
                "canonical PTIRs in one group should have the same number of junctions"
            );

            for (&(left, right), range) in junctions
                .iter()
                .zip(self.canonical_junction_range.iter_mut())
            {
                range.0 = range.0.min(left);
                range.1 = range.1.max(left);
                range.2 = range.2.min(right);
                range.3 = range.3.max(right);
            }
        }

        // select repr_junction and repr_left and repr_right based on policy
        // select representive junction
        for junction_idx in 0..junction_count {
            let positions: Vec<(u32, u32)> = self
                .all_canonical_ptir_list
                .iter()
                .map(|(_, _, _, junc)| junc[junction_idx])
                .collect();

            let (repr, used_policy) = select_pair(&positions, args.splice_policy)?;
            if matches!(args.splice_policy, MergePolicy::Major)
                && matches!(used_policy, MergePolicy::Outer)
            {
                self.used_repr_junction_policy = MergePolicy::Outer;
            }
            self.repr_junction.push(repr);
        }

        // select terminals

        let ((repr_tss, repr_tes), (tss_policy, tes_policy)) = select_repr_terminals(
            &self.all_canonical_ptir_list,
            self.strand,
            args.tss_policy,
            args.tes_policy,
        )?;
        self.set_repr_from_terminals(repr_tss, repr_tes);
        self.set_used_repr_terminal_policies(tss_policy, tes_policy);
        self.repr_loaded = true;
        Ok(())
    }

    pub fn add_non_canonical_ptir(
        &mut self,
        ptir: &PTIR,
        scluster_idx: usize,
    ) -> Result<(), MergeError> {
        self.no_all_canonical_ptir_counts += 1;
        let junction = ptir
            .junction_vec
            .clone()
            .ok_or(MergeError::NoJunctionFound)?;
        self.no_all_canonical_ptir_list
            .push((scluster_idx, ptir.start, ptir.end, junction));
        Ok(())
    }

    pub fn profile_non_canonical_ptirs(&mut self, args: &MergeArgs) -> Result<(), MergeError> {
        self.repr_junction.clear();
        self.repr_left = 0;
        self.repr_right = 0;
        self.used_repr_junction_policy = args.splice_policy;
        self.used_repr_left_policy = args.tss_policy;
        self.used_repr_right_policy = args.tes_policy;
        self.used_repr_mono_policy = args.mono_policy;

        let Some((_, _, _, first_junctions)) = self.no_all_canonical_ptir_list.first() else {
            return Err(MergeError::NoJunctionFound);
        };

        let junction_count = first_junctions.len();

        for (_, _, _, junctions) in self.no_all_canonical_ptir_list.iter().skip(1) {
            debug_assert_eq!(
                junctions.len(),
                junction_count,
                "non-canonical PTIRs in one group should have the same number of junctions"
            );
        }

        for junction_idx in 0..junction_count {
            let positions: Vec<(u32, u32)> = self
                .no_all_canonical_ptir_list
                .iter()
                .map(|(_, _, _, junc)| junc[junction_idx])
                .collect();

            let (repr, used_policy) = select_pair(&positions, args.splice_policy)?;
            if matches!(args.splice_policy, MergePolicy::Major)
                && matches!(used_policy, MergePolicy::Outer)
            {
                self.used_repr_junction_policy = MergePolicy::Outer;
            }
            self.repr_junction.push(repr);
        }

        let ((repr_tss, repr_tes), (tss_policy, tes_policy)) = select_repr_terminals(
            &self.no_all_canonical_ptir_list,
            self.strand,
            args.tss_policy,
            args.tes_policy,
        )?;
        self.set_repr_from_terminals(repr_tss, repr_tes);
        self.set_used_repr_terminal_policies(tss_policy, tes_policy);
        self.repr_loaded = true;
        Ok(())
    }

    pub fn add_mono_ptir(&mut self, ptir: &PTIR, scluster_idx: usize) {
        self.all_canonical_ptir_counts += 1;
        self.all_canonical_ptir_list
            .push((scluster_idx, ptir.start, ptir.end, Vec::new()));
    }

    pub fn profile_mono_ptirs(&mut self, args: &MergeArgs) -> Result<(), MergeError> {
        self.canonical_junction_range.clear();
        self.repr_junction.clear();
        self.repr_left = 0;
        self.repr_right = 0;
        self.used_repr_junction_policy = args.splice_policy;
        self.used_repr_left_policy = args.tss_policy;
        self.used_repr_right_policy = args.tes_policy;
        self.used_repr_mono_policy = args.mono_policy;

        if self.all_canonical_ptir_list.is_empty() {
            return Err(MergeError::SelectReprFailed);
        }

        let positions: Vec<(u32, u32)> = self
            .all_canonical_ptir_list
            .iter()
            .map(|(_, left, right, _)| (*left, *right))
            .collect();

        let ((repr_left, repr_right), mono_policy) = select_pair(&positions, args.mono_policy)?;

        self.repr_left = repr_left;
        self.repr_right = repr_right;
        self.used_repr_mono_policy = mono_policy;
        self.repr_loaded = true;
        Ok(())
    }

    pub fn repr_junction(&self) -> &Vec<(u32, u32)> {
        &self.repr_junction
    }

    pub fn nc_count(&self) -> u32 {
        self.no_all_canonical_ptir_counts
    }

    pub fn ca_count(&self) -> u32 {
        self.all_canonical_ptir_counts
    }

    fn exons_from_repr(&self) -> Vec<(u32, u32)> {
        if self.repr_junction.is_empty() {
            return vec![(self.repr_left, self.repr_right)];
        }

        let mut exons = Vec::with_capacity(self.repr_junction.len() + 1);
        exons.push((self.repr_left, self.repr_junction[0].0));
        for junction_window in self.repr_junction.windows(2) {
            exons.push((junction_window[0].1, junction_window[1].0));
        }
        exons.push((self.repr_junction.last().unwrap().1, self.repr_right));
        exons
    }

    pub fn write_gtf_block(
        &self,
        chrom_name: &str,
        super_cluster: &[PTIR],
        bufwriter: &mut dyn Write,
    ) -> Result<(), MergeError> {
        debug_assert!(
            self.repr_loaded,
            "GroupedPTIR repr must be profiled before GTF export"
        );
        if !self.repr_loaded {
            return Ok(());
        }

        let gene_id = format!("ISOMG_{}", self.gene_id);
        let tx_id = format!("ISOMT_{}", self.tx_id);
        let exons = self.exons_from_repr();
        let tx_type_members = self
            .all_canonical_ptir_list
            .iter()
            .chain(self.no_all_canonical_ptir_list.iter())
            .map(|(super_idx, _, _, _)| &super_cluster[*super_idx])
            .collect::<Vec<_>>();

        let mut source_txs = tx_type_members;
        source_txs.sort_by(|left, right| {
            (
                left.source_file_id,
                left.start,
                left.end,
                left.source_txid.as_str(),
            )
                .cmp(&(
                    right.source_file_id,
                    right.start,
                    right.end,
                    right.source_txid.as_str(),
                ))
        });

        let source_attr = source_txs
            .into_iter()
            .map(|ptir| {
                let (donor_diff, acceptor_diff) = junction_diff_sums(
                    ptir.junctions().unwrap_or(&[]),
                    &self.repr_junction,
                    self.strand,
                );
                format!(
                    "S{}:{}:{}:{}:{}:{}:{}",
                    ptir.source_file_id + 1,
                    ptir.source_txid,
                    ptir.start,
                    ptir.end,
                    tx_type_label(&ptir.tx_type),
                    donor_diff,
                    acceptor_diff
                )
            })
            .collect::<Vec<_>>()
            .join(",");

        let strand = strand_char(self.strand);

        write!(bufwriter, "{chrom_name}\tisomatch\ttranscript\t")?;
        write!(
            bufwriter,
            "{}\t{}\t.\t{}\t.\t",
            self.repr_left, self.repr_right, strand
        )?;
        bufwriter.write_all(b"gene_id \"")?;
        bufwriter.write_all(gene_id.as_bytes())?;
        bufwriter.write_all(b"\"; transcript_id \"")?;
        bufwriter.write_all(tx_id.as_bytes())?;

        bufwriter.write_all(b"\"; ISOM_COUNT \"")?;
        write!(
            bufwriter,
            "{}",
            self.all_canonical_ptir_counts + self.no_all_canonical_ptir_counts
        )?;
        bufwriter.write_all(b"\"; ISOM_SRC \"")?;
        bufwriter.write_all(source_attr.as_bytes())?;

        let isom_policy = if self.n_exon == 1 {
            format!(
                "MONO:NA:NA:NA:{}",
                merge_policy_label(self.used_repr_mono_policy)
            )
        } else {
            format!(
                "MULT:{}:{}:{}:NA",
                merge_policy_label(self.used_repr_junction_policy),
                merge_policy_label(self.used_repr_left_policy),
                merge_policy_label(self.used_repr_right_policy)
            )
        };

        bufwriter.write_all(b"\"; ISOM_POLICY \"")?;
        bufwriter.write_all(isom_policy.as_bytes())?;
        bufwriter.write_all(b"\";\n")?;

        for (idx, (start, end)) in exons.iter().enumerate() {
            write!(bufwriter, "{chrom_name}\tisomatch\texon\t")?;
            write!(bufwriter, "{start}\t{end}\t.\t{strand}\t.\t")?;
            bufwriter.write_all(b"gene_id \"")?;
            bufwriter.write_all(gene_id.as_bytes())?;
            bufwriter.write_all(b"\"; transcript_id \"")?;
            bufwriter.write_all(tx_id.as_bytes())?;
            bufwriter.write_all(b"\"; exon_number \"")?;
            write!(bufwriter, "{}", idx + 1)?;
            bufwriter.write_all(b"\";\n")?;
        }

        Ok(())
    }

    pub fn update_ids(&mut self, gene_id: u32, tx_id: u32) {
        self.gene_id = gene_id;
        self.tx_id = tx_id;
    }
}

impl std::fmt::Display for GroupedPTIR {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            ">>> strand:{}, n_exon:{}, canonical_ptirs:{}, non_canonical_ptirs:{}, repr_start:{}, repr_end:{}, repr_junctions: {:?}
            \ncanonical_ptirs: {:?},
            \nnon canonical ptirs: {:?}
            ",
            self.strand,
            self.n_exon,
            self.all_canonical_ptir_counts,
            self.no_all_canonical_ptir_counts,
            self.repr_left,
            self.repr_right,
            self.repr_junction,
            self.all_canonical_ptir_list,
            self.no_all_canonical_ptir_list
        )
    }
}

fn majority_vote_unique_pair(positions: &[(u32, u32)]) -> Option<(u32, u32)> {
    let mut counts: FxHashMap<(u32, u32), u32> = FxHashMap::default();
    let mut best_count = 0;

    for &position in positions {
        let count = counts.entry(position).or_insert(0);
        *count += 1;
        best_count = best_count.max(*count);
    }

    let mut winners = counts
        .into_iter()
        .filter_map(|(position, count)| (count == best_count).then_some(position));

    let winner = winners.next()?;
    if winners.next().is_some() {
        None
    } else {
        Some(winner)
    }
}

fn outer_pair(positions: &[(u32, u32)]) -> Result<(u32, u32), MergeError> {
    let left = positions
        .iter()
        .map(|(left, _)| *left)
        .min()
        .ok_or(MergeError::SelectReprFailed)?;
    let right = positions
        .iter()
        .map(|(_, right)| *right)
        .max()
        .ok_or(MergeError::SelectReprFailed)?;
    Ok((left, right))
}

fn inner_pair(positions: &[(u32, u32)]) -> Result<(u32, u32), MergeError> {
    let left = positions
        .iter()
        .map(|(left, _)| *left)
        .max()
        .ok_or(MergeError::SelectReprFailed)?;
    let right = positions
        .iter()
        .map(|(_, right)| *right)
        .min()
        .ok_or(MergeError::SelectReprFailed)?;
    Ok((left, right))
}

fn majority_vote_unique_scalar(positions: &[u32]) -> Option<u32> {
    let mut counts: FxHashMap<u32, u32> = FxHashMap::default();
    let mut best_count = 0;

    for &position in positions {
        let count = counts.entry(position).or_insert(0);
        *count += 1;
        best_count = best_count.max(*count);
    }

    let mut winners = counts
        .into_iter()
        .filter_map(|(position, count)| (count == best_count).then_some(position));

    let winner = winners.next()?;
    if winners.next().is_some() {
        None
    } else {
        Some(winner)
    }
}

fn select_pair(
    positions: &[(u32, u32)],
    policy: MergePolicy,
) -> Result<((u32, u32), MergePolicy), MergeError> {
    let out = match policy {
        MergePolicy::Outer => (outer_pair(positions)?, MergePolicy::Outer),
        MergePolicy::Inner => (inner_pair(positions)?, MergePolicy::Inner),
        MergePolicy::Major => match majority_vote_unique_pair(positions) {
            Some(pair) => (pair, MergePolicy::Major),
            None => (outer_pair(positions)?, MergePolicy::Outer),
        },
    };
    Ok(out)
}

fn tss_is_left_boundary(strand: ISOMSTRAND) -> bool {
    match strand {
        ISOMSTRAND::Plus => true,
        ISOMSTRAND::Minus => false,
        ISOMSTRAND::Unknown => true,
    }
}

fn tes_is_left_boundary(strand: ISOMSTRAND) -> bool {
    !tss_is_left_boundary(strand)
}

fn boundaries_to_terminals(left: u32, right: u32, strand: ISOMSTRAND) -> (u32, u32) {
    match strand {
        ISOMSTRAND::Plus => (left, right),
        ISOMSTRAND::Minus => (right, left),
        ISOMSTRAND::Unknown => (left, right),
    }
}

fn terminals_to_boundaries(tss: u32, tes: u32, strand: ISOMSTRAND) -> (u32, u32) {
    match strand {
        ISOMSTRAND::Plus => (tss, tes),
        ISOMSTRAND::Minus => (tes, tss),
        ISOMSTRAND::Unknown => (tss, tes),
    }
}

fn collect_tss_positions(entries: &[GroupedPTIREntry], strand: ISOMSTRAND) -> Vec<u32> {
    entries
        .iter()
        .map(|(_, left, right, _)| boundaries_to_terminals(*left, *right, strand).0)
        .collect()
}

fn collect_tes_positions(entries: &[GroupedPTIREntry], strand: ISOMSTRAND) -> Vec<u32> {
    entries
        .iter()
        .map(|(_, left, right, _)| boundaries_to_terminals(*left, *right, strand).1)
        .collect()
}

fn select_repr_terminals(
    entries: &[GroupedPTIREntry],
    strand: ISOMSTRAND,
    tss_policy: MergePolicy,
    tes_policy: MergePolicy,
) -> Result<((u32, u32), (MergePolicy, MergePolicy)), MergeError> {
    let tss_positions = collect_tss_positions(entries, strand);
    let tes_positions = collect_tes_positions(entries, strand);
    let (repr_tss, used_tss_policy) =
        select_terminal(&tss_positions, tss_policy, tss_is_left_boundary(strand))?;
    let (repr_tes, used_tes_policy) =
        select_terminal(&tes_positions, tes_policy, tes_is_left_boundary(strand))?;
    Ok(((repr_tss, repr_tes), (used_tss_policy, used_tes_policy)))
}

fn select_terminal(
    positions: &[u32],
    policy: MergePolicy,
    is_left_boundary: bool,
) -> Result<(u32, MergePolicy), MergeError> {
    let out = match policy {
        MergePolicy::Outer => {
            if is_left_boundary {
                (
                    *positions.iter().min().ok_or(MergeError::SelectReprFailed)?,
                    MergePolicy::Outer,
                )
            } else {
                (
                    *positions.iter().max().ok_or(MergeError::SelectReprFailed)?,
                    MergePolicy::Outer,
                )
            }
        }
        MergePolicy::Inner => {
            if is_left_boundary {
                (
                    *positions.iter().max().ok_or(MergeError::SelectReprFailed)?,
                    MergePolicy::Inner,
                )
            } else {
                (
                    *positions.iter().min().ok_or(MergeError::SelectReprFailed)?,
                    MergePolicy::Inner,
                )
            }
        }
        MergePolicy::Major => match majority_vote_unique_scalar(positions) {
            Some(pos) => (pos, MergePolicy::Major),
            None => {
                if is_left_boundary {
                    (
                        *positions.iter().min().ok_or(MergeError::SelectReprFailed)?,
                        MergePolicy::Outer,
                    )
                } else {
                    (
                        *positions.iter().max().ok_or(MergeError::SelectReprFailed)?,
                        MergePolicy::Outer,
                    )
                }
            }
        },
    };
    Ok(out)
}

fn junction_diff_sums(curr: &[(u32, u32)], repr: &[(u32, u32)], strand: ISOMSTRAND) -> (u32, u32) {
    if curr.len() != repr.len() {
        return (u32::MAX, u32::MAX);
    }

    let mut donor_sum = 0;
    let mut acceptor_sum = 0;
    for (curr_junc, repr_junc) in curr.iter().zip(repr.iter()) {
        match strand {
            ISOMSTRAND::Plus => {
                donor_sum += curr_junc.0.abs_diff(repr_junc.0);
                acceptor_sum += curr_junc.1.abs_diff(repr_junc.1);
            }
            ISOMSTRAND::Minus => {
                donor_sum += curr_junc.1.abs_diff(repr_junc.1);
                acceptor_sum += curr_junc.0.abs_diff(repr_junc.0);
            }
            ISOMSTRAND::Unknown => {
                donor_sum += curr_junc.0.abs_diff(repr_junc.0);
                acceptor_sum += curr_junc.1.abs_diff(repr_junc.1);
            }
        }
    }
    (donor_sum, acceptor_sum)
}

fn strand_char(strand: ISOMSTRAND) -> char {
    match strand {
        ISOMSTRAND::Plus => '+',
        ISOMSTRAND::Minus => '-',
        ISOMSTRAND::Unknown => '.',
    }
}

fn tx_type_label(tx_type: &TxType) -> &'static str {
    match tx_type {
        TxType::MONO => "MONO",
        TxType::ALLC => "ALL_CA",
        TxType::PRTC => "PRT_CA",
        TxType::NOTC => "NOT_CA",
    }
}

fn merge_policy_label(policy: MergePolicy) -> &'static str {
    match policy {
        MergePolicy::Outer => "outer",
        MergePolicy::Inner => "inner",
        MergePolicy::Major => "major",
    }
}
