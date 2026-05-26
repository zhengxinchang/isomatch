use std::io::Write;

use crate::{
    MergeArgs,
    core::{ptir::PTIR, tx_strand::ISOMSTRAND, tx_type::TxType},
    merge::{
        guide::GuideDb,
        merge_error::MergeError,
        policy::{GuideResolution, MergePolicyArg, MergePolicyUsed},
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
    // used_repr_mono_policy: MergePolicyUsed,
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
            // used_repr_mono_policy: MergePolicyUsed::Major,
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

    pub fn from_canonical_entries(
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
            // used_repr_mono_policy: MergePolicyUsed::Major,
            repr_loaded: false,
        }
    }

    pub fn from_non_canonical_entries(
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
            // used_repr_mono_policy: MergePolicyUsed::Major,
            repr_loaded: false,
        }
    }

    fn set_repr_from_terminals(&mut self, tss: u32, tes: u32) {
        let (repr_left, repr_right) = match self.strand {
            ISOMSTRAND::Plus => (tss, tes),
            ISOMSTRAND::Minus => (tes, tss),
            ISOMSTRAND::Unknown => (tss, tes),
        };
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
        // self.used_repr_mono_policy = MergePolicyUsed::from_arg_policy(&args.mono_policy);

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

            let (repr, used_policy) = select_splice_pair(&positions, args.splice_policy)?;
            if matches!(args.splice_policy, MergePolicyArg::Major)
                && matches!(used_policy, MergePolicyArg::Longer)
            {
                self.used_repr_junction_policy = MergePolicyUsed::Longer;
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
        // self.used_repr_mono_policy = MergePolicyUsed::from_arg_policy(&args.mono_policy);

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

            let (repr, used_policy) = select_splice_pair(&positions, args.splice_policy)?;
            if matches!(args.splice_policy, MergePolicyArg::Major)
                && matches!(used_policy, MergePolicyArg::Longer)
            {
                self.used_repr_junction_policy = MergePolicyUsed::Longer;
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

    pub fn profile_mono_ptirs(
        &mut self,
        chrom: &str,
        args: &MergeArgs,
        guide_tss: &Option<GuideDb>,
        guide_tes: &Option<GuideDb>,
    ) -> Result<(), MergeError> {
        self.canonical_junction_range.clear();
        self.repr_junction.clear();
        self.repr_left = 0;
        self.repr_right = 0;
        self.used_repr_junction_policy = MergePolicyUsed::from_arg_policy(&args.splice_policy);
        self.used_repr_left_policy = MergePolicyUsed::from_arg_policy(&args.tss_policy);
        self.used_repr_right_policy = MergePolicyUsed::from_arg_policy(&args.tes_policy);
        // self.used_repr_mono_policy = MergePolicyUsed::from_arg_policy(&args.mono_policy);

        if self.all_canonical_ptir_list.is_empty() {
            return Err(MergeError::SelectReprFailed);
        }

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

    pub fn repr_junction(&self) -> &Vec<(u32, u32)> {
        &self.repr_junction
    }

    pub fn nc_count(&self) -> u32 {
        self.no_all_canonical_ptir_counts
    }

    pub fn ca_count(&self) -> u32 {
        self.all_canonical_ptir_counts
    }

    pub fn total_count(&self) -> u32 {
        self.all_canonical_ptir_counts + self.no_all_canonical_ptir_counts
    }

    pub fn used_tss_policy(&self) -> MergePolicyUsed {
        if tss_is_left_boundary(&self.strand) {
            self.used_repr_left_policy
        } else {
            self.used_repr_right_policy
        }
    }

    pub fn used_tes_policy(&self) -> MergePolicyUsed {
        if tss_is_left_boundary(&self.strand) {
            self.used_repr_right_policy
        } else {
            self.used_repr_left_policy
        }
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
        gtf_bufwriter: &mut dyn Write,
        track_bufwriter: &mut dyn Write,
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

        let src_records: Vec<SrcRecord> = source_txs
            .iter()
            .map(|ptir| {
                let (donor_diff, acceptor_diff) = junction_diff_sums(
                    ptir.junctions().unwrap_or(&[]),
                    &self.repr_junction,
                    self.strand,
                );
                let exons_diff =
                    junction_exon_diffs(ptir.junctions().unwrap_or(&[]), &self.repr_junction)
                        .expect("PTIR must have same junctions as representative");

                let exon_diff_str = if exons_diff.is_empty() {
                    "no_diff".to_string()
                } else {
                    exons_diff
                        .into_iter()
                        .map(|a| format!("({},{},{})", a.0, a.1, a.2))
                        .collect::<Vec<_>>()
                        .join(",")
                };

                let gtf_str = format!(
                    "S{}:{}:{}:{}:{}:{}:{}:{}",
                    ptir.source_file_id + 1,
                    ptir.source_txid,
                    ptir.start,
                    ptir.end,
                    ptir.tx_type,
                    donor_diff,
                    acceptor_diff,
                    exon_diff_str
                );

                SrcRecord {
                    gtf_str,
                    ptir,
                    donor_diff,
                    acceptor_diff,
                    exon_diff_str,
                }
            })
            .collect();

        let source_attr = src_records
            .iter()
            .map(|r| r.gtf_str.as_str())
            .collect::<Vec<_>>()
            .join("|");

        let strand = char::from(self.strand);

        write!(gtf_bufwriter, "{chrom_name}\tisomatch\ttranscript\t")?;
        write!(
            gtf_bufwriter,
            "{}\t{}\t.\t{}\t.\t",
            self.repr_left, self.repr_right, strand
        )?;
        gtf_bufwriter.write_all(b"gene_id \"")?;
        gtf_bufwriter.write_all(gene_id.as_bytes())?;
        gtf_bufwriter.write_all(b"\"; transcript_id \"")?;
        gtf_bufwriter.write_all(tx_id.as_bytes())?;

        let isom_exons = self.n_exon.to_string();
        gtf_bufwriter.write_all(b"\"; ISOM_EXONS \"")?;
        gtf_bufwriter.write_all(isom_exons.as_bytes())?;

        gtf_bufwriter.write_all(b"\"; ISOM_COUNT \"")?;
        write!(
            gtf_bufwriter,
            "{}",
            self.all_canonical_ptir_counts + self.no_all_canonical_ptir_counts
        )?;

        gtf_bufwriter.write_all(b"\"; ISOM_SRC \"")?;
        gtf_bufwriter.write_all(source_attr.as_bytes())?;

        let isom_policy = format!(
            "{}:{}:{}",
            if self.n_exon == 1 {
                "NA".to_string()
            } else {
                self.used_repr_junction_policy.to_string()
            },
            self.used_repr_left_policy,
            self.used_repr_right_policy,
        );

        gtf_bufwriter.write_all(b"\"; ISOM_REPR_POLICY \"")?;
        gtf_bufwriter.write_all(isom_policy.as_bytes())?;
        gtf_bufwriter.write_all(b"\";\n")?;

        for (idx, (start, end)) in exons.iter().enumerate() {
            write!(gtf_bufwriter, "{chrom_name}\tisomatch\texon\t")?;
            write!(gtf_bufwriter, "{start}\t{end}\t.\t{strand}\t.\t")?;
            gtf_bufwriter.write_all(b"gene_id \"")?;
            gtf_bufwriter.write_all(gene_id.as_bytes())?;
            gtf_bufwriter.write_all(b"\"; transcript_id \"")?;
            gtf_bufwriter.write_all(tx_id.as_bytes())?;
            gtf_bufwriter.write_all(b"\"; exon_number \"")?;
            write!(gtf_bufwriter, "{}", idx + 1)?;
            gtf_bufwriter.write_all(b"\";\n")?;
        }

        let total_src_count = self.all_canonical_ptir_counts + self.no_all_canonical_ptir_counts;
        let junction_policy = if self.n_exon == 1 {
            "NA".to_string()
        } else {
            self.used_repr_junction_policy.to_string()
        };
        let tss_policy = self.used_tss_policy().to_string();
        let tes_policy = self.used_tes_policy().to_string();
        for r in &src_records {
            writeln!(
                track_bufwriter,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                tx_id,
                gene_id,
                self.repr_left,
                self.repr_right,
                strand,
                self.n_exon,
                junction_policy,
                tss_policy,
                tes_policy,
                total_src_count,
                r.ptir.source_txid,
                r.ptir.source_geneid,
                r.donor_diff,
                r.acceptor_diff,
                r.exon_diff_str,
            )?;
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

/// temp output struct for src item
struct SrcRecord<'a> {
    gtf_str: String,
    ptir: &'a PTIR,
    donor_diff: u32,
    acceptor_diff: u32,
    exon_diff_str: String,
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

fn select_splice_pair(
    positions: &[(u32, u32)],
    policy: MergePolicyArg,
) -> Result<((u32, u32), MergePolicyArg), MergeError> {
    let out = match policy {
        // For splice junctions, a shorter intron yields longer flanking exons.
        MergePolicyArg::Longer => (inner_pair(positions)?, MergePolicyArg::Longer),
        MergePolicyArg::Shorter => (outer_pair(positions)?, MergePolicyArg::Shorter),
        MergePolicyArg::Major => match majority_vote_unique_pair(positions) {
            Some(pair) => (pair, MergePolicyArg::Major),
            None => (inner_pair(positions)?, MergePolicyArg::Longer),
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

fn boundaries_to_terminals(left: u32, right: u32, strand: ISOMSTRAND) -> (u32, u32) {
    match strand {
        ISOMSTRAND::Plus => (left, right),
        ISOMSTRAND::Minus => (right, left),
        ISOMSTRAND::Unknown => (left, right),
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

fn majority_vote_unique_position(positions: &[u32]) -> Option<u32> {
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

fn select_terminal_by_policy(
    positions: &[u32],
    is_left_boundary: bool,
    policy: MergePolicyArg,
) -> Result<(u32, MergePolicyUsed), MergeError> {
    let choose_longer = |positions: &[u32]| -> Result<u32, MergeError> {
        if is_left_boundary {
            positions
                .iter()
                .copied()
                .min()
                .ok_or(MergeError::SelectReprFailed)
        } else {
            positions
                .iter()
                .copied()
                .max()
                .ok_or(MergeError::SelectReprFailed)
        }
    };

    let choose_shorter = |positions: &[u32]| -> Result<u32, MergeError> {
        if is_left_boundary {
            positions
                .iter()
                .copied()
                .max()
                .ok_or(MergeError::SelectReprFailed)
        } else {
            positions
                .iter()
                .copied()
                .min()
                .ok_or(MergeError::SelectReprFailed)
        }
    };

    match policy {
        MergePolicyArg::Longer => Ok((choose_longer(positions)?, MergePolicyUsed::Longer)),
        MergePolicyArg::Shorter => Ok((choose_shorter(positions)?, MergePolicyUsed::Shorter)),
        MergePolicyArg::Major => match majority_vote_unique_position(positions) {
            Some(position) => Ok((position, MergePolicyUsed::Major)),
            None => Ok((choose_longer(positions)?, MergePolicyUsed::Longer)),
        },
    }
}

fn select_terminal(
    chrom: &str,
    positions: &[u32],
    strand: &ISOMSTRAND,
    is_left_boundary: bool,
    policy: MergePolicyArg,
    guide: &Option<GuideDb>,
    guide_flank: u32,
) -> Result<(u32, MergePolicyUsed), MergeError> {
    // if no guide file provided, fallback to non guide policy
    let Some(guide) = guide.as_ref() else {
        return select_terminal_by_policy(positions, is_left_boundary, policy);
    };

    // for each tss or tes, check how many of guide regions ovlp with it
    // and record the max hits
    let mut max_hits = 0usize;
    let mut hits = Vec::with_capacity(positions.len());
    for &position in positions {
        let hit_count = guide
            .query_overlaps_with_flank(chrom, strand, position, guide_flank)
            .len();
        max_hits = max_hits.max(hit_count);
        hits.push(hit_count);
    }

    // if no guide region overlapped with any tss/tes position, fallback to non-guide policy
    if max_hits == 0 {
        return select_terminal_by_policy(positions, is_left_boundary, policy);
    }

    // for pick up the max ovlpped tss/tes
    let guide_candidates: Vec<u32> = positions
        .iter()
        .zip(hits.iter())
        .filter_map(|(&position, &hit_count)| (hit_count == max_hits).then_some(position))
        .collect();

    // guard assert, make sure at lest one element exists
    if guide_candidates.is_empty() {
        return Err(MergeError::SelectReprFailed);
    }

    // major vote if guide_candidates
    // if all candidates point to the same position, guide is the sole determinant
    // e.g. guide_candidates = [100,100,100], each one have 2 guide region support
    if guide_candidates.iter().all(|&c| c == guide_candidates[0]) {
        return Ok((
            guide_candidates[0],
            MergePolicyUsed::Guide(GuideResolution::Definitive),
        ));
    }

    // else ==> find the most frequent positions and return
    if let Some(position) = majority_vote_unique_position(&guide_candidates) {
        return Ok((position, MergePolicyUsed::Guide(GuideResolution::Majority)));
    }

    // at this step, the guide candidates have same guide ovlps and they have same
    // recurrency in the list.
    // the only choice is choose the longer one.
    Ok((
        select_terminal_by_policy(&guide_candidates, is_left_boundary, MergePolicyArg::Longer)?.0,
        MergePolicyUsed::Guide(GuideResolution::Longer),
    ))
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
    let (repr_tss, used_tss) = select_terminal(
        chrom,
        &tss_pos,
        strand,
        tss_is_left_boundary(strand),
        tss_policy,
        guide_tss,
        guide_tss_flank,
    )?;
    let (repr_tes, used_tes) = select_terminal(
        chrom,
        &tes_pos,
        strand,
        !tss_is_left_boundary(strand),
        tes_policy,
        guide_tes,
        guide_tes_flank,
    )?;

    Ok(((repr_tss, repr_tes), (used_tss, used_tes)))
}

/// calculate the sum of junction difference between current transcript and repr.
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
        i32,   // left diff bp
        i32,   // right diff bp
    )>,
    MergeError,
> {
    if curr.len() != repr.len() {
        return Err(MergeError::NoSameJunction);
    }
    let mut exon_diffs = Vec::new();
    for (exon_idx, (curr_junc, repr_junc)) in curr.iter().zip(repr.iter()).enumerate() {
        let left_diff = repr_junc.0 as i32 - curr_junc.0 as i32;
        let right_diff = repr_junc.1 as i32 - curr_junc.1 as i32;
        if left_diff > 0 || right_diff > 0 {
            exon_diffs.push((exon_idx + 1, left_diff, right_diff))
        }
    }
    Ok(exon_diffs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splice_longer_prefers_shorter_intron_span() {
        let positions = vec![(100, 200), (110, 190), (105, 195)];
        let (repr, used_policy) = select_splice_pair(&positions, MergePolicyArg::Longer).unwrap();
        assert_eq!(repr, (110, 190));
        assert!(matches!(used_policy, MergePolicyArg::Longer));
    }
}
