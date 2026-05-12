use std::io::Write;

use crate::{
    MergeArgs,
    core::{ptir::PTIR, status::TxType, tx_strand::ISOMSTRAND},
    merge::{
        guide::GuideDb,
        merge_error::MergeError,
        policy::{MergePolicyArg, MergePolicyUsed},
    },
};
use rustc_hash::FxHashMap;

#[derive(Clone, Debug)]
pub struct GroupedPTIREntry {
    pub super_idx: usize,
    pub left: u32,
    pub right: u32,
    pub junctions: Vec<(u32, u32)>,
    pub tx_type: TxType,
}
impl GroupedPTIREntry {
    pub fn tss(&self, strand: &ISOMSTRAND) -> u32 {
        match strand {
            ISOMSTRAND::Plus | ISOMSTRAND::Unknown => self.left,
            ISOMSTRAND::Minus => self.right,
        }
    }

    pub fn tes(&self, strand: &ISOMSTRAND) -> u32 {
        match strand {
            ISOMSTRAND::Plus | ISOMSTRAND::Unknown => self.right,
            ISOMSTRAND::Minus => self.left,
        }
    }
}

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
    used_repr_junction_policy: MergePolicyUsed,
    used_repr_left_policy: MergePolicyUsed,
    used_repr_right_policy: MergePolicyUsed,
    used_repr_mono_policy: MergePolicyUsed,
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
            used_repr_junction_policy: MergePolicyUsed::Major,
            used_repr_left_policy: MergePolicyUsed::Major,
            used_repr_right_policy: MergePolicyUsed::Major,
            used_repr_mono_policy: MergePolicyUsed::Major,
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
            used_repr_junction_policy: MergePolicyUsed::Major,
            used_repr_left_policy: MergePolicyUsed::Major,
            used_repr_right_policy: MergePolicyUsed::Major,
            used_repr_mono_policy: MergePolicyUsed::Major,
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
            used_repr_junction_policy: MergePolicyUsed::Major,
            used_repr_left_policy: MergePolicyUsed::Major,
            used_repr_right_policy: MergePolicyUsed::Major,
            used_repr_mono_policy: MergePolicyUsed::Major,
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
        tss_policy: MergePolicyUsed,
        tes_policy: MergePolicyUsed,
    ) {
        if tss_is_left_boundary(&self.strand) {
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
        self.all_canonical_ptir_list.push(GroupedPTIREntry {
            super_idx: scluster_idx,
            left: ptir.start,
            right: ptir.end,
            junctions: junction,
            tx_type: ptir.tx_type,
        });
        Ok(())
    }

    pub fn profile_canonical_ptirs(
        &mut self,
        chrom: &str,
        args: &MergeArgs,
        guide_tss: &Option<GuideDb>,
        guide_tes: &Option<GuideDb>,
    ) -> Result<(), MergeError> {
        // generate canonical_junction_range based on canonical_ptir_list
        self.canonical_junction_range.clear();
        self.repr_junction.clear();
        self.repr_left = 0;
        self.repr_right = 0;
        self.used_repr_junction_policy = MergePolicyUsed::from_arg_policy(&args.splice_policy);
        self.used_repr_left_policy = MergePolicyUsed::from_arg_policy(&args.tss_policy);
        self.used_repr_right_policy = MergePolicyUsed::from_arg_policy(&args.tes_policy);
        self.used_repr_mono_policy = MergePolicyUsed::from_arg_policy(&args.mono_policy);

        let Some(first_entry) = self.all_canonical_ptir_list.first() else {
            return Err(MergeError::NoJunctionFound);
        };

        let junction_count = first_entry.junctions.len();
        self.canonical_junction_range = first_entry
            .junctions
            .iter()
            .map(|&(left, right)| (left, left, right, right))
            .collect();

        for entry in self.all_canonical_ptir_list.iter().skip(1) {
            debug_assert_eq!(
                entry.junctions.len(),
                self.canonical_junction_range.len(),
                "canonical PTIRs in one group should have the same number of junctions"
            );

            for (&(left, right), range) in entry
                .junctions
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
                .map(|entry| entry.junctions[junction_idx])
                .collect();

            let (repr, used_policy) = select_pair(&positions, args.splice_policy)?;
            if matches!(args.splice_policy, MergePolicyArg::Major)
                && matches!(used_policy, MergePolicyArg::Outer)
            {
                self.used_repr_junction_policy = MergePolicyUsed::Outer;
            }
            self.repr_junction.push(repr);
        }

        // select terminals

        let ((repr_tss, repr_tes), (tss_policy, tes_policy)) = select_repr_terminals(
            chrom,
            &self.all_canonical_ptir_list,
            &self.strand,
            args.tss_policy,
            args.tes_policy,
            guide_tss,
            guide_tes,
            args.guide_tss_flank,
            args.guide_tes_flank,
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
        self.no_all_canonical_ptir_list.push(GroupedPTIREntry {
            super_idx: scluster_idx,
            left: ptir.start,
            right: ptir.end,
            junctions: junction,
            tx_type: ptir.tx_type,
        });
        Ok(())
    }

    pub fn profile_non_canonical_ptirs(
        &mut self,
        chrom: &str,
        args: &MergeArgs,
        guide_tss: &Option<GuideDb>,
        guide_tes: &Option<GuideDb>,
    ) -> Result<(), MergeError> {
        self.repr_junction.clear();
        self.repr_left = 0;
        self.repr_right = 0;
        self.used_repr_junction_policy = MergePolicyUsed::from_arg_policy(&args.splice_policy);
        self.used_repr_left_policy = MergePolicyUsed::from_arg_policy(&args.tss_policy);
        self.used_repr_right_policy = MergePolicyUsed::from_arg_policy(&args.tes_policy);
        self.used_repr_mono_policy = MergePolicyUsed::from_arg_policy(&args.mono_policy);

        let Some(first_entry) = self.no_all_canonical_ptir_list.first() else {
            return Err(MergeError::NoJunctionFound);
        };

        let junction_count = first_entry.junctions.len();

        for entry in self.no_all_canonical_ptir_list.iter().skip(1) {
            debug_assert_eq!(
                entry.junctions.len(),
                junction_count,
                "non-canonical PTIRs in one group should have the same number of junctions"
            );
        }

        for junction_idx in 0..junction_count {
            let positions: Vec<(u32, u32)> = self
                .no_all_canonical_ptir_list
                .iter()
                .map(|entry| entry.junctions[junction_idx])
                .collect();

            let (repr, used_policy) = select_pair(&positions, args.splice_policy)?;
            if matches!(args.splice_policy, MergePolicyArg::Major)
                && matches!(used_policy, MergePolicyArg::Outer)
            {
                self.used_repr_junction_policy = MergePolicyUsed::Outer;
            }
            self.repr_junction.push(repr);
        }

        let ((repr_tss, repr_tes), (tss_policy, tes_policy)) = select_repr_terminals(
            chrom,
            &self.no_all_canonical_ptir_list,
            &self.strand,
            args.tss_policy,
            args.tes_policy,
            guide_tss,
            guide_tes,
            args.guide_tss_flank,
            args.guide_tes_flank,
        )?;
        self.set_repr_from_terminals(repr_tss, repr_tes);
        self.set_used_repr_terminal_policies(tss_policy, tes_policy);
        self.repr_loaded = true;
        Ok(())
    }

    pub fn add_mono_ptir(&mut self, ptir: &PTIR, scluster_idx: usize) {
        self.all_canonical_ptir_counts += 1;
        self.all_canonical_ptir_list.push(GroupedPTIREntry {
            super_idx: scluster_idx,
            left: ptir.start,
            right: ptir.end,
            junctions: Vec::new(),
            tx_type: ptir.tx_type,
        });
    }

    pub fn profile_mono_ptirs(&mut self, args: &MergeArgs) -> Result<(), MergeError> {
        self.canonical_junction_range.clear();
        self.repr_junction.clear();
        self.repr_left = 0;
        self.repr_right = 0;
        self.used_repr_junction_policy = MergePolicyUsed::from_arg_policy(&args.splice_policy);
        self.used_repr_left_policy = MergePolicyUsed::from_arg_policy(&args.tss_policy);
        self.used_repr_right_policy = MergePolicyUsed::from_arg_policy(&args.tes_policy);
        self.used_repr_mono_policy = MergePolicyUsed::from_arg_policy(&args.mono_policy);

        if self.all_canonical_ptir_list.is_empty() {
            return Err(MergeError::SelectReprFailed);
        }

        let positions: Vec<(u32, u32)> = self
            .all_canonical_ptir_list
            .iter()
            .map(|entry| (entry.left, entry.right))
            .collect();

        let ((repr_left, repr_right), mono_policy) = select_pair(&positions, args.mono_policy)?;

        self.repr_left = repr_left;
        self.repr_right = repr_right;
        self.used_repr_mono_policy = MergePolicyUsed::from_arg_policy(&mono_policy);
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
            .map(|entry| &super_cluster[entry.super_idx])
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
                let exons_diff =
                    junction_exon_diffs(ptir.junctions().unwrap_or(&[]), &self.repr_junction)
                        .expect("PTIR must have same junctions as representative");

                let exon_diff_string = if exons_diff.len() > 0 {
                    exons_diff
                        .into_iter()
                        .map(|a| format!("({},{},{})", a.0, a.1, a.2))
                        .collect::<Vec<_>>()
                        .join(",")
                } else {
                    "no_diff".to_string()
                };

                format!(
                    "S{}:{}:{}:{}:{}:{}:{}:{}",
                    ptir.source_file_id + 1,
                    ptir.source_txid,
                    ptir.start,
                    ptir.end,
                    ptir.tx_type,
                    donor_diff,
                    acceptor_diff,
                    exon_diff_string
                )
            })
            .collect::<Vec<_>>()
            .join("|");

        let strand = char::from(self.strand);

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

        let isom_exons = if self.n_exon == 1 { "MONO" } else { "MULTI" };
        bufwriter.write_all(b"\"; ISOM_EXONS \"")?;
        bufwriter.write_all(isom_exons.as_bytes())?;

        bufwriter.write_all(b"\"; ISOM_COUNT \"")?;
        write!(
            bufwriter,
            "{}",
            self.all_canonical_ptir_counts + self.no_all_canonical_ptir_counts
        )?;

        bufwriter.write_all(b"\"; ISOM_SRC \"")?;
        bufwriter.write_all(source_attr.as_bytes())?;

        let isom_policy = if self.n_exon == 1 {
            format!("NA:NA:NA:{}", self.used_repr_mono_policy)
        } else {
            format!(
                "{}:{}:{}:NA",
                self.used_repr_junction_policy,
                self.used_repr_left_policy,
                self.used_repr_right_policy
            )
        };

        bufwriter.write_all(b"\"; ISOM_REPR_POLICY \"")?;
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

fn select_pair(
    positions: &[(u32, u32)],
    policy: MergePolicyArg,
) -> Result<((u32, u32), MergePolicyArg), MergeError> {
    let out = match policy {
        MergePolicyArg::Outer => (outer_pair(positions)?, MergePolicyArg::Outer),
        MergePolicyArg::Inner => (inner_pair(positions)?, MergePolicyArg::Inner),
        MergePolicyArg::Major => match majority_vote_unique_pair(positions) {
            Some(pair) => (pair, MergePolicyArg::Major),
            None => (outer_pair(positions)?, MergePolicyArg::Outer),
        },
    };
    Ok(out)
}

fn tss_is_left_boundary(strand: &ISOMSTRAND) -> bool {
    match strand {
        ISOMSTRAND::Plus => true,
        ISOMSTRAND::Minus => false,
        ISOMSTRAND::Unknown => true,
    }
}

fn tes_is_left_boundary(strand: &ISOMSTRAND) -> bool {
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

fn collect_tss_positions(entries: &[GroupedPTIREntry], strand: &ISOMSTRAND) -> Vec<u32> {
    entries
        .iter()
        .map(|entry| boundaries_to_terminals(entry.left, entry.right, *strand).0)
        .collect()
}

fn collect_tes_positions(entries: &[GroupedPTIREntry], strand: &ISOMSTRAND) -> Vec<u32> {
    entries
        .iter()
        .map(|entry| boundaries_to_terminals(entry.left, entry.right, *strand).1)
        .collect()
}

fn select_repr_terminals(
    chrom: &str,
    entries: &[GroupedPTIREntry],
    strand: &ISOMSTRAND,
    tss_policy: MergePolicyArg,
    tes_policy: MergePolicyArg,
    guide_tss: &Option<GuideDb>,
    guide_tes: &Option<GuideDb>,
    guide_tss_flank: u32,
    guide_tes_flank: u32,
) -> Result<((u32, u32), (MergePolicyUsed, MergePolicyUsed)), MergeError> {
    let tss_pos = collect_tss_positions(entries, strand);
    let tes_pos = collect_tes_positions(entries, strand);

    // Step 1: compute guide hit counts for every entry up front
    let (tss_hits, tes_hits): (Vec<usize>, Vec<usize>) = (0..entries.len())
        .map(|i| {
            let t = guide_tss.as_ref().map_or(0, |g| {
                g.query_overlaps_with_flank(chrom, strand, tss_pos[i], guide_tss_flank)
                    .len()
            });
            let e = guide_tes.as_ref().map_or(0, |g| {
                g.query_overlaps_with_flank(chrom, strand, tes_pos[i], guide_tes_flank)
                    .len()
            });
            (t, e)
        })
        .unzip();

    // Step 2: classify entries by guide support
    let full_guided: Vec<usize> = (0..entries.len())
        .filter(|&i| {
            guide_tss.is_some() && guide_tes.is_some() && tss_hits[i] > 0 && tes_hits[i] > 0
        })
        .collect();

    let partial_guided: Vec<usize> = (0..entries.len())
        .filter(|&i| tss_hits[i] > 0 || tes_hits[i] > 0)
        .collect();

    // Step 3: select representative in priority order

    // Case A: at least one entry has both TSS and TES guide support
    //         → use position policy only as a tiebreaker among them
    if !full_guided.is_empty() {
        let idx = select_repr_by_policy(
            &tss_pos,
            &tes_pos,
            &full_guided,
            strand,
            tss_policy,
            tes_policy,
        )?
        .0;
        return Ok((
            (tss_pos[idx], tes_pos[idx]),
            (MergePolicyUsed::Guide, MergePolicyUsed::Guide),
        ));
    }

    // Case B: some entries have partial guide support (TSS-only or TES-only)
    //         → pick by total hit score, tie-break by transcript length, then input order
    if !partial_guided.is_empty() {
        let idx = select_repr_by_guide_score(entries, &tss_hits, &tes_hits, &partial_guided);
        let used_tss = if tss_hits[idx] > 0 {
            MergePolicyUsed::Guide
        } else {
            MergePolicyUsed::from_arg_policy(&tss_policy)
        };
        let used_tes = if tes_hits[idx] > 0 {
            MergePolicyUsed::Guide
        } else {
            MergePolicyUsed::from_arg_policy(&tes_policy)
        };
        return Ok(((tss_pos[idx], tes_pos[idx]), (used_tss, used_tes)));
    }

    // Case C: no guide support → fall back to position-based policy on all entries
    let all: Vec<usize> = (0..entries.len()).collect();
    let (idx, used_tss, used_tes) =
        select_repr_by_policy(&tss_pos, &tes_pos, &all, strand, tss_policy, tes_policy)?;
    Ok(((tss_pos[idx], tes_pos[idx]), (used_tss, used_tes)))
}

fn select_repr_by_guide_score(
    entries: &[GroupedPTIREntry],
    tss_hits: &[usize],
    tes_hits: &[usize],
    candidates: &[usize],
) -> usize {
    let mut best = candidates[0];
    for &idx in candidates.iter().skip(1) {
        let best_score = tss_hits[best] + tes_hits[best];
        let curr_score = tss_hits[idx] + tes_hits[idx];
        if curr_score > best_score {
            best = idx;
            continue;
        }
        if curr_score == best_score {
            let best_len = entries[best].right - entries[best].left + 1;
            let curr_len = entries[idx].right - entries[idx].left + 1;
            if curr_len > best_len {
                best = idx;
            }
        }
        // tied on score and length → keep earliest (current best unchanged)
    }
    best
}

fn select_repr_by_policy(
    tss_pos: &[u32],
    tes_pos: &[u32],
    candidates: &[usize],
    strand: &ISOMSTRAND,
    tss_policy: MergePolicyArg,
    tes_policy: MergePolicyArg,
) -> Result<(usize, MergePolicyUsed, MergePolicyUsed), MergeError> {
    if candidates.is_empty() {
        return Err(MergeError::SelectReprFailed);
    }

    let mut tss_counts: FxHashMap<u32, usize> = FxHashMap::default();
    let mut tes_counts: FxHashMap<u32, usize> = FxHashMap::default();
    for &idx in candidates {
        *tss_counts.entry(tss_pos[idx]).or_insert(0) += 1;
        *tes_counts.entry(tes_pos[idx]).or_insert(0) += 1;
    }

    let max_tss_count = tss_counts.values().copied().max().unwrap_or(0);
    let max_tes_count = tes_counts.values().copied().max().unwrap_or(0);
    let unique_tss_major = tss_counts.values().filter(|&&c| c == max_tss_count).count() == 1;
    let unique_tes_major = tes_counts.values().filter(|&&c| c == max_tes_count).count() == 1;

    let used_tss_policy = match tss_policy {
        MergePolicyArg::Outer => MergePolicyUsed::Outer,
        MergePolicyArg::Inner => MergePolicyUsed::Inner,
        MergePolicyArg::Major => {
            if unique_tss_major {
                MergePolicyUsed::Major
            } else {
                MergePolicyUsed::Outer
            }
        }
    };
    let used_tes_policy = match tes_policy {
        MergePolicyArg::Outer => MergePolicyUsed::Outer,
        MergePolicyArg::Inner => MergePolicyUsed::Inner,
        MergePolicyArg::Major => {
            if unique_tes_major {
                MergePolicyUsed::Major
            } else {
                MergePolicyUsed::Outer
            }
        }
    };

    let mut best = candidates[0];
    for &idx in candidates.iter().skip(1) {
        let best_tss = tss_pos[best];
        let curr_tss = tss_pos[idx];
        let best_tss_count = *tss_counts.get(&best_tss).unwrap_or(&0);
        let curr_tss_count = *tss_counts.get(&curr_tss).unwrap_or(&0);

        // curr_better_tss and best_better_tss are mutually exclusive, so the
        // negation guards on the other side are redundant and omitted.
        let curr_better_tss = match tss_policy {
            MergePolicyArg::Outer => {
                if tss_is_left_boundary(strand) {
                    curr_tss < best_tss
                } else {
                    curr_tss > best_tss
                }
            }
            MergePolicyArg::Inner => {
                if tss_is_left_boundary(strand) {
                    curr_tss > best_tss
                } else {
                    curr_tss < best_tss
                }
            }
            MergePolicyArg::Major => {
                if curr_tss_count != best_tss_count {
                    curr_tss_count > best_tss_count
                } else if tss_is_left_boundary(strand) {
                    curr_tss < best_tss
                } else {
                    curr_tss > best_tss
                }
            }
        };
        let best_better_tss = match tss_policy {
            MergePolicyArg::Outer => {
                if tss_is_left_boundary(strand) {
                    best_tss < curr_tss
                } else {
                    best_tss > curr_tss
                }
            }
            MergePolicyArg::Inner => {
                if tss_is_left_boundary(strand) {
                    best_tss > curr_tss
                } else {
                    best_tss < curr_tss
                }
            }
            MergePolicyArg::Major => {
                if best_tss_count != curr_tss_count {
                    best_tss_count > curr_tss_count
                } else if tss_is_left_boundary(strand) {
                    best_tss < curr_tss
                } else {
                    best_tss > curr_tss
                }
            }
        };

        if curr_better_tss {
            best = idx;
            continue;
        }
        if best_better_tss {
            continue;
        }

        // TSS tied → compare TES; if also tied, keep earliest (best unchanged)
        let best_tes = tes_pos[best];
        let curr_tes = tes_pos[idx];
        let best_tes_count = *tes_counts.get(&best_tes).unwrap_or(&0);
        let curr_tes_count = *tes_counts.get(&curr_tes).unwrap_or(&0);

        let curr_better_tes = match tes_policy {
            MergePolicyArg::Outer => {
                if tes_is_left_boundary(strand) {
                    curr_tes < best_tes
                } else {
                    curr_tes > best_tes
                }
            }
            MergePolicyArg::Inner => {
                if tes_is_left_boundary(strand) {
                    curr_tes > best_tes
                } else {
                    curr_tes < best_tes
                }
            }
            MergePolicyArg::Major => {
                if curr_tes_count != best_tes_count {
                    curr_tes_count > best_tes_count
                } else if tes_is_left_boundary(strand) {
                    curr_tes < best_tes
                } else {
                    curr_tes > best_tes
                }
            }
        };

        if curr_better_tes {
            best = idx;
        }
    }

    Ok((best, used_tss_policy, used_tes_policy))
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

fn junction_exon_diffs(
    curr: &[(u32, u32)],
    repr: &[(u32, u32)],
) -> Result<
    Vec<(
        usize, // exon number
        u32,   // left diff bp
        u32,   // right diff bp
    )>,
    MergeError,
> {
    if curr.len() != repr.len() {
        return Err(MergeError::NoSameJunction);
    }
    let mut exon_diffs = Vec::new();
    for (exon_idx, (curr_junc, repr_junc)) in curr.iter().zip(repr.iter()).enumerate() {
        let left_diff = repr_junc.0 - curr_junc.0;
        let right_diff = repr_junc.1 - curr_junc.1;
        if left_diff > 0 || right_diff > 0 {
            exon_diffs.push((exon_idx + 1, left_diff, right_diff))
        }
    }
    Ok(exon_diffs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::merge::guide::GuideBEDType;
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn unique_temp_path(prefix: &str, suffix: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("{prefix}-{nanos}.{suffix}"))
    }

    fn load_guide_db(contents: &str, guide_type: GuideBEDType) -> GuideDb {
        let path = unique_temp_path("isomatch-grouped-ptirs-guide", "bed");
        fs::write(&path, contents).unwrap();
        let chrmap: Option<&PathBuf> = None;
        let db = GuideDb::from_bed_path(&path, guide_type, &chrmap).unwrap();
        let _ = fs::remove_file(&path);
        db
    }

    fn entry(left: u32, right: u32) -> GroupedPTIREntry {
        GroupedPTIREntry {
            super_idx: 0,
            left,
            right,
            junctions: Vec::new(),
            tx_type: TxType::MONO,
        }
    }

    #[test]
    fn select_repr_terminals_prefers_transcript_supported_on_both_sides() {
        let entries = vec![entry(100, 200), entry(120, 220)];
        let guide_tss = Some(load_guide_db(
            concat!(
                "chromosome\tstart\tend\tID\tscore\tstrand\n",
                "chr1\t99\t100\ttss0\t1\t+\n",
            ),
            GuideBEDType::Tss,
        ));
        let guide_tes = Some(load_guide_db(
            concat!(
                "chromosome\tstart\tend\tID\tscore\tstrand\n",
                "chr1\t199\t200\ttes0\t1\t+\n",
                "chr1\t219\t220\ttes1\t1\t+\n",
                "chr1\t218\t221\ttes2\t1\t+\n",
            ),
            GuideBEDType::Tes,
        ));

        let selected = select_repr_terminals(
            "chr1",
            &entries,
            &ISOMSTRAND::Plus,
            MergePolicyArg::Outer,
            MergePolicyArg::Outer,
            &guide_tss,
            &guide_tes,
            0,
            0,
        )
        .unwrap();

        assert_eq!(selected.0, (100, 200));
        assert!(matches!(
            selected.1,
            (MergePolicyUsed::Guide, MergePolicyUsed::Guide)
        ));
    }

    #[test]
    fn select_repr_terminals_uses_policy_within_multi_guided_candidates() {
        let entries = vec![entry(100, 200), entry(110, 190)];
        let guide_tss = Some(load_guide_db(
            concat!(
                "chromosome\tstart\tend\tID\tscore\tstrand\n",
                "chr1\t99\t100\ttss0\t1\t+\n",
                "chr1\t109\t110\ttss1\t1\t+\n",
            ),
            GuideBEDType::Tss,
        ));
        let guide_tes = Some(load_guide_db(
            concat!(
                "chromosome\tstart\tend\tID\tscore\tstrand\n",
                "chr1\t199\t200\ttes0\t1\t+\n",
                "chr1\t189\t190\ttes1\t1\t+\n",
            ),
            GuideBEDType::Tes,
        ));

        let selected = select_repr_terminals(
            "chr1",
            &entries,
            &ISOMSTRAND::Plus,
            MergePolicyArg::Outer,
            MergePolicyArg::Outer,
            &guide_tss,
            &guide_tes,
            0,
            0,
        )
        .unwrap();

        assert_eq!(selected.0, (100, 200));
        assert!(matches!(
            selected.1,
            (MergePolicyUsed::Guide, MergePolicyUsed::Guide)
        ));
    }

    #[test]
    fn select_repr_terminals_scores_partial_support_then_breaks_ties_by_length() {
        let entries = vec![entry(100, 180), entry(120, 230)];
        let guide_tss = Some(load_guide_db(
            concat!(
                "chromosome\tstart\tend\tID\tscore\tstrand\n",
                "chr1\t99\t100\ttss0\t1\t+\n",
                "chr1\t119\t120\ttss1\t1\t+\n",
            ),
            GuideBEDType::Tss,
        ));
        let guide_tes: Option<GuideDb> = None;

        let selected = select_repr_terminals(
            "chr1",
            &entries,
            &ISOMSTRAND::Plus,
            MergePolicyArg::Outer,
            MergePolicyArg::Inner,
            &guide_tss,
            &guide_tes,
            0,
            0,
        )
        .unwrap();

        assert_eq!(selected.0, (120, 230));
        assert!(matches!(
            selected.1,
            (MergePolicyUsed::Guide, MergePolicyUsed::Inner)
        ));
    }

    #[test]
    fn select_repr_terminals_without_guide_support_falls_back_to_policy_on_real_transcript() {
        let entries = vec![entry(100, 220), entry(110, 200)];
        let guide_tss: Option<GuideDb> = None;
        let guide_tes: Option<GuideDb> = None;

        let selected = select_repr_terminals(
            "chr1",
            &entries,
            &ISOMSTRAND::Plus,
            MergePolicyArg::Outer,
            MergePolicyArg::Inner,
            &guide_tss,
            &guide_tes,
            0,
            0,
        )
        .unwrap();

        assert_eq!(selected.0, (100, 220));
        assert!(matches!(
            selected.1,
            (MergePolicyUsed::Outer, MergePolicyUsed::Inner)
        ));
    }

    #[test]
    fn select_repr_terminals_uses_earliest_transcript_for_full_ties() {
        let entries = vec![entry(100, 200), entry(120, 220)];
        let guide_tss = Some(load_guide_db(
            concat!(
                "chromosome\tstart\tend\tID\tscore\tstrand\n",
                "chr1\t99\t100\ttss0\t1\t+\n",
                "chr1\t119\t120\ttss1\t1\t+\n",
            ),
            GuideBEDType::Tss,
        ));
        let guide_tes: Option<GuideDb> = None;

        let selected = select_repr_terminals(
            "chr1",
            &entries,
            &ISOMSTRAND::Plus,
            MergePolicyArg::Outer,
            MergePolicyArg::Outer,
            &guide_tss,
            &guide_tes,
            0,
            0,
        )
        .unwrap();

        assert_eq!(selected.0, (100, 200));
        assert!(matches!(
            selected.1,
            (MergePolicyUsed::Guide, MergePolicyUsed::Outer)
        ));
    }
}
