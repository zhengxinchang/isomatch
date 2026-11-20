use std::collections::{BTreeMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
pub mod bench_error;
use crate::merge::policy::{MergePolicyArg, TerminalRefineMode};
use anyhow::{Context, Result as AnyResult, bail};
use bench_error::BenchError;
use clap::ValueEnum;
use flate2::bufread::GzDecoder;
use flate2::{Compression, write::GzEncoder};
use libgtf::gtf::{Strand, attribute, parse_record};
use log::{error, info, warn};
use serde::Serialize;

use crate::constants::ISOM_GTF_SCHEMA;
use crate::index::IndexStats;
use crate::merge::run_merge;
use crate::{
    BenchArgs, MergeArgs,
    tools::revert::{parse_junction_diffs, shift_coordinate, write_gtf_record},
    traits::ArgValidate,
    utils::{greetings2, parse_header_value, require_file, save_json_block},
};

/// `--preset`: what a match means and which records are evaluated.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, ValueEnum)]
#[serde(rename_all = "lowercase")]
pub enum BenchPreset {
    /// Original records; junctions and TSS/TES are compared.
    Default,
    /// Each input collapsed by exact junction chain (mono: exact ends); only junctions compared.
    Junction,
}

const DEFAULT_END_WOB: u32 = 50;

impl BenchArgs {
    /// TSS/TES thresholds of the default preset; None for junction (ends are not compared).
    fn end_wobs(&self) -> Option<(u32, u32)> {
        (self.preset == BenchPreset::Default).then(|| {
            (
                self.tss_wob.unwrap_or(DEFAULT_END_WOB),
                self.tes_wob.unwrap_or(DEFAULT_END_WOB),
            )
        })
    }
}

impl ArgValidate for BenchArgs {
    fn validate(&self) {
        let mut error_msg = String::new();
        let mut has_error = false;

        require_file(
            "Reference GTF",
            &self.ref_gtf,
            &mut error_msg,
            &mut has_error,
        );
        require_file(
            "Reference FASTA",
            &self.ref_fa,
            &mut error_msg,
            &mut has_error,
        );

        let mut ref_fai = self.ref_fa.clone();
        ref_fai.add_extension("fai");
        require_file(
            "Reference FASTA index",
            &ref_fai,
            &mut error_msg,
            &mut has_error,
        );

        if !(self.mono_ovlp.is_finite() && self.mono_ovlp > 0.0 && self.mono_ovlp <= 1.0) {
            error_msg.push_str("\n--mono-ovlp must be in (0, 1].");
            has_error = true;
        }

        if self.preset == BenchPreset::Junction
            && (self.tss_wob.is_some() || self.tes_wob.is_some())
        {
            error_msg.push_str("\n--tss-wob/--tes-wob only apply to --preset default.");
            has_error = true;
        }

        // Existing-directory rules (empty, marker, --overwrite) are handled by prepare_outdir.
        if self.outdir.exists() && !self.outdir.is_dir() {
            error_msg.push_str("\nOutput path is not a directory.");
            has_error = true;
        }

        if has_error {
            error!("Error validating arguments: {}", error_msg);
            std::process::exit(1);
        }
    }
}

pub fn run_bench(args: BenchArgs) -> AnyResult<()> {
    greetings2(&args);
    args.validate();
    let mut ref_fai = args.ref_fa.clone();
    ref_fai.add_extension("fai");
    let mut protected: Vec<PathBuf> = args.inputs.clone();
    protected.extend([args.ref_gtf.clone(), args.ref_fa.clone(), ref_fai]);
    protected.push(std::env::current_dir()?);
    let mut forbidden = vec![PathBuf::from("/")];
    forbidden.extend(std::env::var_os("HOME").map(PathBuf::from));
    prepare_outdir(&args.outdir, args.overwrite, &protected, &forbidden)?;

    let mut param_path = args.outdir.clone();
    param_path.push("params.json");
    save_json_block(&param_path, &args)?;

    // Merge sample order: S1 = reference (sample_idx 0), S2.. = queries (sample_idx 1..).
    // User-facing outputs label them REF and Q1.. instead (see `sample_label`).
    let mut raw_inputs = args.inputs.clone();
    raw_inputs.insert(0, args.ref_gtf.clone());

    // Evaluated GTFs: the inputs themselves, or (junction) each input collapsed on its own by
    // merge with zero wobble and mono_ovlp 1.0, i.e. exact junction chain / exact mono ends.
    let eval_paths = match args.preset {
        BenchPreset::Default => raw_inputs.clone(),
        BenchPreset::Junction => raw_inputs
            .iter()
            .enumerate()
            .map(|(idx, raw)| {
                let label = sample_label(idx);
                let dir = args.outdir.join("collapsed").join(&label);
                std::fs::create_dir_all(&dir)?;
                let prefix = dir.join(&label);
                info!(
                    "Collapsing {label} by exact junction chain: {}",
                    raw.display()
                );
                run_merge(bench_merge_args(
                    &args,
                    vec![raw.clone()],
                    prefix.clone(),
                    0,
                    1.0,
                ))?;
                let mut gtf = prefix;
                gtf.add_extension("merged.gtf.gz");
                Ok(gtf)
            })
            .collect::<AnyResult<Vec<_>>>()?,
    };

    // Merge output is an intermediate grouping step, kept under initial_group/.
    let group_dir = args.outdir.join("initial_group");
    std::fs::create_dir_all(&group_dir)?;
    // run_merge checks inputs and builds missing/outdated indexes (collapsed GTFs included).
    info!("Merging transcripts from {} input GTFs", eval_paths.len());
    run_merge(bench_merge_args(
        &args,
        eval_paths.clone(),
        group_dir.join("bench"),
        args.wob,
        args.mono_ovlp,
    ))?;

    let inputs = raw_inputs
        .into_iter()
        .zip(eval_paths)
        .map(|(raw_path, eval_path)| {
            let raw_stats = load_index_stats(&raw_path, args.skip_missing_ref_chr)?;
            let collapsed = match args.preset {
                BenchPreset::Default => None,
                BenchPreset::Junction => Some(CollapsedInput {
                    gtf_stats: load_index_stats(&eval_path, args.skip_missing_ref_chr)?,
                    path: eval_path,
                }),
            };
            Ok(BenchInput {
                raw_path,
                raw_stats,
                collapsed,
            })
        })
        .collect::<AnyResult<Vec<_>>>()?;
    let totals: Vec<u64> = inputs
        .iter()
        .map(|input| input.eval_stats().transcript_count)
        .collect();

    info!("Matching query transcripts to reference transcripts");
    let mut recovered = vec![0u64; inputs.len()];
    let mut matches = vec![
        SampleMatches {
            multi: vec![0; args.wob as usize + 1],
            mono: 0,
        };
        inputs.len()
    ];
    // Per query sample status GTFs at w = 0; index 0 (reference) has no writers.
    let mut status_writers = (1..inputs.len())
        .map(|sample_idx| open_status_writers(&args.outdir, sample_idx, args.preset, &inputs))
        .collect::<AnyResult<Vec<_>>>()?;
    let merged_path = group_dir.join("bench.merged.gtf.gz");
    for group in TxReader::new(&merged_path)? {
        let group = group?;
        for tx in &group.tx_structures {
            let Some(count) = recovered.get_mut(tx.sample_idx) else {
                bail!(
                    "Sample S{} in {} is out of range",
                    tx.sample_idx + 1,
                    merged_path.display()
                );
            };
            *count += 1;
        }
        let w0_pairs = match_group(&group, &args, &mut matches);
        write_group_status(&group, &w0_pairs, &mut status_writers)?;
    }
    for writers in status_writers {
        for writer in writers {
            writer.into_inner()?.finish()?;
        }
    }

    for ((input, &seen), &total) in inputs.iter().zip(&recovered).zip(&totals) {
        if seen != total {
            bail!(
                "Recovered {seen} transcripts for {} from merged GTF, but its index reports {total}",
                input.eval_path().display()
            );
        }
    }

    let report = build_report(args.preset, inputs, &matches)?;
    let report_path = args.outdir.join("results_all.json");
    save_json_block(&report_path, &report)?;
    let summary_path = args.outdir.join("results_summary.json");
    save_json_block(&summary_path, &build_summary(&report, &args))?;
    let tsv_path = args.outdir.join("results_summary.tsv");
    write_summary_tsv(&report, &tsv_path)?;
    for q in &report.queries {
        let m = &q.results[0].all;
        info!(
            "{} (w=0): TP={} FP={} FN={} precision={} recall={} F1={}",
            q.input.sample_id,
            m.tp_query,
            m.fp,
            m.fn_,
            fmt_ratio(m.precision),
            fmt_ratio(m.recall),
            fmt_ratio(m.f1)
        );
    }
    info!("Bench report saved to: {}", report_path.display());
    info!("Bench summary saved to: {}", summary_path.display());
    info!("Bench summary table saved to: {}", tsv_path.display());
    info!("Finished!");
    Ok(())
}

fn fmt_ratio(value: Option<f64>) -> String {
    value.map_or("NA".to_string(), |v| v.to_string())
}

/// results_summary.tsv: one row per query and w, metrics for all records.
fn write_summary_tsv(report: &BenchReport, path: &Path) -> AnyResult<()> {
    let mut out =
        String::from("sample_id\tinput_path\tw\ttp_query\ttp_ref\tfp\tfn\tprecision\trecall\tf1\n");
    for q in &report.queries {
        for r in &q.results {
            let m = &r.all;
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                q.input.sample_id,
                q.input.input_path.display(),
                r.w,
                m.tp_query,
                m.tp_ref,
                m.fp,
                m.fn_,
                fmt_ratio(m.precision),
                fmt_ratio(m.recall),
                fmt_ratio(m.f1)
            ));
        }
    }
    std::fs::write(path, out).with_context(|| format!("Failed to write {}", path.display()))
}

/// MergeArgs for both merges in bench: the per-input collapse (wob 0, mono_ovlp 1.0) and the
/// grouping merge (`args.wob`, `args.mono_ovlp`). All six splice wobbles are set explicitly
/// even though ignore_canonical makes merge ignore the *_nc ones; terminal ends never split.
fn bench_merge_args(
    args: &BenchArgs,
    inputs: Vec<PathBuf>,
    out: PathBuf,
    wob: u32,
    mono_ovlp: f64,
) -> MergeArgs {
    // Unused with terminal_refine None; kept equal to the bench thresholds for the params file.
    let (tss_wob, tes_wob) = args
        .end_wobs()
        .unwrap_or((DEFAULT_END_WOB, DEFAULT_END_WOB));
    MergeArgs {
        inputs,
        out,
        chop: false,
        ref_fa: args.ref_fa.clone(),
        skip_missing_ref_chr: args.skip_missing_ref_chr,
        ignore_canonical: true,
        wob_d: wob,
        wob_a: wob,
        wob_u: wob,
        terminal_refine: TerminalRefineMode::None,
        tss_wob,
        tes_wob,
        wob_d_nc: wob,
        wob_a_nc: wob,
        wob_u_nc: wob,
        terminal_refine_nc: TerminalRefineMode::None,
        tss_wob_nc: tss_wob,
        tes_wob_nc: tes_wob,
        mono_ovlp,
        guide_tss: None,
        guide_tes: None,
        guide_tss_flank: 10,
        guide_tes_flank: 10,
        chrmap: None,
        splice_policy: MergePolicyArg::Major,
        tss_policy: MergePolicyArg::Major,
        tes_policy: MergePolicyArg::Major,
        gene_assign_min_overlap: 80,
        unassigned_group_min_overlap: 20,
        quiet: true,
    }
}

/// One bench input: the command-line GTF and, for the junction preset, its collapsed GTF.
struct BenchInput {
    raw_path: PathBuf,
    raw_stats: IndexStats,
    collapsed: Option<CollapsedInput>,
}

impl BenchInput {
    /// GTF whose records are evaluated.
    fn eval_path(&self) -> &Path {
        self.collapsed.as_ref().map_or(&self.raw_path, |c| &c.path)
    }

    /// Index summary of the evaluated GTF; its counts are the denominators.
    fn eval_stats(&self) -> &IndexStats {
        self.collapsed
            .as_ref()
            .map_or(&self.raw_stats, |c| &c.gtf_stats)
    }
}

#[derive(Debug, Serialize)]
struct CollapsedInput {
    path: PathBuf,
    gtf_stats: IndexStats,
}

/// Matched pairs of one query sample: multi-exon per w = 0..=wob; mono-exon does not use w.
#[derive(Debug, Clone, Default)]
struct SampleMatches {
    multi: Vec<u64>,
    mono: u64,
}

/// Marker written into every bench output directory; only marked directories are cleared.
const BENCH_MARKER: &str = ".ISOMATCH_BENCH";

/// Create or reuse the output directory and write the marker:
/// - missing or empty: used as is;
/// - non-empty without --overwrite: error;
/// - non-empty with --overwrite: contents are deleted only if the marker is present, no
///   `protected` path (inputs, reference, working directory) lies inside it and it is not one
///   of `forbidden` (root, home). Symlinks inside are removed, not followed.
fn prepare_outdir(
    outdir: &Path,
    overwrite: bool,
    protected: &[PathBuf],
    forbidden: &[PathBuf],
) -> AnyResult<()> {
    let is_empty = match std::fs::read_dir(outdir) {
        Ok(mut entries) => entries.next().is_none(),
        Err(err) if err.kind() == std::io::ErrorKind::NotFound => true,
        Err(err) => {
            return Err(err)
                .with_context(|| format!("Can not read output directory {}", outdir.display()));
        }
    };
    if !is_empty {
        if !overwrite {
            bail!(
                "Output directory {} is not empty, use --overwrite to replace a previous bench output",
                outdir.display()
            );
        }
        if !outdir.join(BENCH_MARKER).is_file() {
            bail!(
                "Refusing to clear output directory {}: it has no {BENCH_MARKER} marker, so it was not created by isomatch bench. Choose another output directory or empty it manually.",
                outdir.display()
            );
        }
        let outdir_canon = outdir.canonicalize()?;
        // Missing paths are reported later by input validation and can not be inside outdir.
        let canon = |path: &PathBuf| path.canonicalize().ok();
        if let Some(path) = protected
            .iter()
            .find(|path| canon(path).is_some_and(|c| c.starts_with(&outdir_canon)))
        {
            bail!(
                "Refusing to clear output directory {}: it contains {}",
                outdir.display(),
                path.display()
            );
        }
        if forbidden
            .iter()
            .any(|path| canon(path) == Some(outdir_canon.clone()))
        {
            bail!(
                "Refusing to clear output directory {}: it is the root or home directory",
                outdir.display()
            );
        }
        warn!(
            "--overwrite: deleting existing contents of {}",
            outdir.display()
        );
        for entry in std::fs::read_dir(outdir)? {
            let path = entry?.path();
            if path.symlink_metadata()?.is_dir() {
                std::fs::remove_dir_all(&path)
            } else {
                std::fs::remove_file(&path)
            }
            .with_context(|| format!("Can not delete {}", path.display()))?;
        }
    }
    std::fs::create_dir_all(outdir)?;
    std::fs::write(
        outdir.join(BENCH_MARKER),
        format!(
            "{} {} bench output directory\n",
            env!("CARGO_PKG_NAME"),
            env!("CARGO_PKG_VERSION")
        ),
    )?;
    Ok(())
}

/// Read the index summary written by `isomatch index` (`<gtf>.index_info.json`).
fn load_index_stats(gtf: &Path, skip_missing_ref_chr: bool) -> AnyResult<IndexStats> {
    let mut info_path = gtf.to_path_buf();
    info_path.add_extension("index_info.json");
    let stats: IndexStats = File::open(&info_path)
        .map_err(anyhow::Error::from)
        .and_then(|f| Ok(serde_json::from_reader(BufReader::new(f))?))
        .with_context(|| {
            format!(
                "Can not read index summary {}, please rerun isomatch index",
                info_path.display()
            )
        })?;
    if stats.missing_seqid_cnt > 0 {
        if !skip_missing_ref_chr {
            bail!(
                "Index of {} skipped transcripts on seqid(s) absent from the reference FASTA: {}. Rerun with --skip-missing-ref-chr to accept this.",
                gtf.display(),
                stats.missing_seqids.join(",")
            );
        }
        warn!(
            "{}: {} transcript(s) on missing reference seqid(s) are excluded",
            gtf.display(),
            stats.skipped_transcript_cnt
        );
    }
    Ok(stats)
}

/// Greedy one-to-one matching per query sample inside one merged group.
/// Adds the number of pairs to `matches[sample_idx]` (multi-exon: for each w = 0..=args.wob;
/// mono-exon: once) and returns the w = 0 pairs as (query, reference) indices into
/// `group.tx_structures`.
fn match_group(
    group: &TxGroup,
    args: &BenchArgs,
    matches: &mut [SampleMatches],
) -> Vec<(usize, usize)> {
    let mut w0_pairs = Vec::new();
    if group.ref_tx_cnt == 0 {
        return w0_pairs;
    }
    let mut refs = Vec::new();
    let mut ref_idxs = Vec::new();
    let mut queries: BTreeMap<usize, (Vec<&TxStructure>, Vec<usize>)> = BTreeMap::new();
    for (idx, tx) in group.tx_structures.iter().enumerate() {
        if tx.is_ref_tx() {
            refs.push(tx);
            ref_idxs.push(idx);
        } else {
            let (qs, q_idxs) = queries.entry(tx.sample_idx).or_default();
            qs.push(tx);
            q_idxs.push(idx);
        }
    }

    for (sample_idx, (qs, q_idxs)) in queries {
        // ponytail: all Q x R pairs in the group are compared, O(Q * R * L); index candidates
        // if large groups become a bottleneck.
        let mut cands = Vec::new();
        for (qi, q) in qs.iter().enumerate() {
            for (ri, r) in refs.iter().enumerate() {
                if let Some(diffs) = candidate_diffs(q, r, args) {
                    cands.push((diffs, qi, ri));
                }
            }
        }
        // Mono-exon candidates have jc_diff_max = 0, so the key reduces to (tss, tes, ids).
        // Junction preset, multi-exon: ends do not take part, so the key is (jc, ids).
        let ignore_ends = args.preset == BenchPreset::Junction && group.n_exons > 1;
        let key = |d: &TxDiffs| {
            if ignore_ends {
                (d.jc_diff_max, 0, 0)
            } else {
                (d.jc_diff_max, d.tss_diff, d.tes_diff)
            }
        };
        cands.sort_by(|(a, aq, ar), (b, bq, br)| {
            (key(a), &qs[*aq].orig_txid, &refs[*ar].orig_txid).cmp(&(
                key(b),
                &qs[*bq].orig_txid,
                &refs[*br].orig_txid,
            ))
        });

        let sample_matches = &mut matches[sample_idx];
        let pairs_w0 = if group.n_exons == 1 {
            let pairs = greedy_match(&cands, qs.len(), refs.len());
            sample_matches.mono += pairs.len() as u64;
            pairs
        } else {
            let mut pairs_w0 = Vec::new();
            for (w, m) in sample_matches.multi.iter_mut().enumerate() {
                // Candidates are sorted by jc_diff_max first, so those within w form a prefix.
                let n_within = cands.partition_point(|(d, _, _)| d.jc_diff_max as usize <= w);
                let pairs = greedy_match(&cands[..n_within], qs.len(), refs.len());
                *m += pairs.len() as u64;
                if w == 0 {
                    pairs_w0 = pairs;
                }
            }
            pairs_w0
        };
        w0_pairs.extend(
            pairs_w0
                .into_iter()
                .map(|(qi, ri)| (q_idxs[qi], ref_idxs[ri])),
        );
    }
    w0_pairs
}

/// Accept candidates in order when both sides are unused; return the accepted (query, ref) pairs.
fn greedy_match(
    cands: &[(TxDiffs, usize, usize)],
    n_query: usize,
    n_ref: usize,
) -> Vec<(usize, usize)> {
    let mut query_used = vec![false; n_query];
    let mut ref_used = vec![false; n_ref];
    let mut pairs = Vec::new();
    for &(_, qi, ri) in cands {
        if !query_used[qi] && !ref_used[ri] {
            query_used[qi] = true;
            ref_used[ri] = true;
            pairs.push((qi, ri));
        }
    }
    pairs
}

/// Status files written under `<outdir>/S#/`, in the order of `StatusWriters`.
const STATUS_FILES: [&str; 4] = ["tp-query", "fp", "tp-ref", "fn"];
type StatusWriters = [BufWriter<GzEncoder<File>>; 4];
const TP_QUERY: usize = 0;
const FP: usize = 1;
const TP_REF: usize = 2;
const FN: usize = 3;

/// Open `<outdir>/Q#/<status>.gtf.gz` for one query sample and write the header of each file.
/// `inputs[0]` is the reference.
fn open_status_writers(
    outdir: &Path,
    sample_idx: usize,
    preset: BenchPreset,
    inputs: &[BenchInput],
) -> AnyResult<StatusWriters> {
    let dir = outdir.join(sample_label(sample_idx));
    std::fs::create_dir_all(&dir)?;
    let open = |name: &str| -> AnyResult<BufWriter<GzEncoder<File>>> {
        let path = dir.join(format!("{name}.gtf.gz"));
        let file =
            File::create(&path).with_context(|| format!("Can not create {}", path.display()))?;
        Ok(BufWriter::new(GzEncoder::new(file, Compression::default())))
    };
    let [tp_query, fp, tp_ref, fn_] = STATUS_FILES.map(open);
    let mut writers = [tp_query?, fp?, tp_ref?, fn_?];
    for (status, writer) in writers.iter_mut().enumerate() {
        write_status_header(writer, status, sample_idx, preset, inputs)?;
    }
    Ok(writers)
}

/// Label shown to users: REF for the reference, Q1.. for queries in command-line order.
/// Merge-internal S# IDs only appear in initial_group/.
fn sample_label(sample_idx: usize) -> String {
    if sample_idx == 0 {
        "REF".to_string()
    } else {
        format!("Q{sample_idx}")
    }
}

/// Header in the merge GTF style; only the STATUS_BENCH tag and the SAMPLE raw_input/role/total
/// fields are bench specific. SAMPLE ids use REF/Q# labels instead of merge S# IDs; `input` is
/// the evaluated GTF and `raw_input` (junction only) the command-line GTF it was collapsed from.
fn write_status_header(
    writer: &mut dyn Write,
    status: usize,
    sample_idx: usize,
    preset: BenchPreset,
    inputs: &[BenchInput],
) -> AnyResult<()> {
    let escape = |value: &str| value.replace('\\', "\\\\").replace('"', "\\\"");
    let sample_id = sample_label(sample_idx);
    let (side, description) = match status {
        TP_QUERY => (
            "query",
            format!("query transcripts of {sample_id} with a one-to-one reference match at w=0"),
        ),
        FP => (
            "query",
            format!("query transcripts of {sample_id} without a one-to-one reference match at w=0"),
        ),
        TP_REF => (
            "reference",
            format!("reference transcripts matched one-to-one by an {sample_id} transcript at w=0"),
        ),
        _ => (
            "reference",
            format!("reference transcripts not matched by any {sample_id} transcript at w=0"),
        ),
    };
    writeln!(
        writer,
        "##ISOM <VERSION> version=\"{}\"; program=\"{}\"; schema=\"{}\"",
        env!("CARGO_PKG_VERSION"),
        env!("CARGO_PKG_NAME"),
        ISOM_GTF_SCHEMA
    )?;
    writeln!(
        writer,
        "##ISOM <STATUS_BENCH> id=\"{}\"; side=\"{side}\"; w=\"0\"; preset=\"{}\"; description=\"{description}\";",
        STATUS_FILES[status],
        match preset {
            BenchPreset::Default => "default",
            BenchPreset::Junction => "junction",
        }
    )?;
    for (idx, role) in [(0, "reference"), (sample_idx, "query")] {
        let input = &inputs[idx];
        let raw_input = match input.collapsed {
            Some(_) => format!(
                " raw_input=\"{}\";",
                escape(&input.raw_path.to_string_lossy())
            ),
            None => String::new(),
        };
        writeln!(
            writer,
            "##ISOM <SAMPLE> id=\"{}\"; input=\"{}\";{raw_input} role=\"{role}\"; total=\"{}\";",
            sample_label(idx),
            escape(&input.eval_path().to_string_lossy()),
            input.eval_stats().transcript_count
        )?;
    }
    writeln!(
        writer,
        "##ISOM <FORMAT> ID=\"ISOM_TX_ID\"; Description=\"merged transcript ID in initial_group/bench.merged.gtf.gz of the group holding this transcript\";"
    )?;
    let command = std::env::args_os()
        .map(|arg| escape(&arg.to_string_lossy()))
        .collect::<Vec<_>>()
        .join(" ");
    writeln!(writer, "##ISOM <COMMAND> cmd={command}")?;
    Ok(())
}

/// Write every record of the group once per query sample: query records of that sample go to
/// tp-query/fp and reference records go to tp-ref/fn, judged by the w = 0 pairs of that sample.
fn write_group_status(
    group: &TxGroup,
    w0_pairs: &[(usize, usize)],
    status_writers: &mut [StatusWriters],
) -> AnyResult<()> {
    // (query sample, record index) for both sides of every pair.
    let mut matched = HashSet::new();
    for &(qi, ri) in w0_pairs {
        let sample_idx = group.tx_structures[qi].sample_idx;
        matched.insert((sample_idx, qi));
        matched.insert((sample_idx, ri));
    }
    for (writer_idx, writers) in status_writers.iter_mut().enumerate() {
        let sample_idx = writer_idx + 1;
        for (idx, tx) in group.tx_structures.iter().enumerate() {
            let is_matched = matched.contains(&(sample_idx, idx));
            let status = match (tx.is_ref_tx(), is_matched) {
                (true, true) => TP_REF,
                (true, false) => FN,
                _ if tx.sample_idx != sample_idx => continue,
                (false, true) => TP_QUERY,
                (false, false) => FP,
            };
            write_tx_gtf(&mut writers[status], group, tx)?;
        }
    }
    Ok(())
}

/// Write one source transcript in the same layout as `isomatch tools revert`.
fn write_tx_gtf(writer: &mut dyn Write, group: &TxGroup, tx: &TxStructure) -> AnyResult<()> {
    let strand = char::from(group.strand).to_string();
    let exons = tx.exons();
    let (start, end) = (exons[0].0, exons[exons.len() - 1].1);
    let record = |writer: &mut dyn Write, feature, start, end, exon_number| {
        write_gtf_record(
            writer,
            &group.chr,
            feature,
            start,
            end,
            &strand,
            &tx.orig_geneid,
            &tx.orig_txid,
            &group.tx_id,
            exon_number,
        )
    };
    record(writer, "transcript", start, end, None)?;
    for (idx, &(exon_start, exon_end)) in exons.iter().enumerate() {
        record(writer, "exon", exon_start, exon_end, Some(idx + 1))?;
    }
    Ok(())
}

/// Return the diffs when a query/reference pair is a candidate under the maximum thresholds.
fn candidate_diffs(
    query: &TxStructure,
    reference: &TxStructure,
    args: &BenchArgs,
) -> Option<TxDiffs> {
    let (q_left, q_right) = (query.tss.min(query.tes), query.tss.max(query.tes));
    let (r_left, r_right) = (
        reference.tss.min(reference.tes),
        reference.tss.max(reference.tes),
    );
    let overlap_left = q_left.max(r_left);
    let overlap_right = q_right.min(r_right);
    if overlap_left > overlap_right {
        return None;
    }

    let diffs = query.compare(reference);
    if let Some((tss_wob, tes_wob)) = args.end_wobs()
        && (diffs.tss_diff > tss_wob || diffs.tes_diff > tes_wob)
    {
        return None;
    }

    if query.junctions.is_empty() {
        let overlap_len = u64::from(overlap_right - overlap_left) + 1;
        let q_len = u64::from(q_right - q_left) + 1;
        let r_len = u64::from(r_right - r_left) + 1;
        if (overlap_len as f64 / q_len as f64) < args.mono_ovlp
            || (overlap_len as f64 / r_len as f64) < args.mono_ovlp
        {
            return None;
        }
    } else if diffs.jc_diff_max > args.wob {
        return None;
    }
    Some(diffs)
}

fn ratio(num: u64, den: u64) -> Option<f64> {
    (den > 0).then(|| num as f64 / den as f64)
}

fn f1_score(precision: Option<f64>, recall: Option<f64>) -> Option<f64> {
    let (p, r) = (precision?, recall?);
    Some(if p + r == 0.0 {
        0.0
    } else {
        2.0 * p * r / (p + r)
    })
}

/// Keep six decimals in the report; `None` (zero denominator) is written as null.
fn round6(value: Option<f64>) -> Option<f64> {
    value.map(|v| (v * 1e6).round() / 1e6)
}

#[derive(Debug, Serialize)]
struct BenchReport {
    preset: BenchPreset,
    reference: InputReport,
    query_count: usize,
    queries: Vec<QueryReport>,
}

/// `total` counts the evaluated records; `gtf_stats` is the command-line GTF's index summary
/// and `collapsed` (junction only, otherwise null) the collapsed GTF and its summary.
#[derive(Debug, Serialize)]
struct InputReport {
    sample_id: String,
    input_path: PathBuf,
    total: u64,
    gtf_stats: IndexStats,
    collapsed: Option<CollapsedInput>,
}

#[derive(Debug, Serialize)]
struct QueryReport {
    #[serde(flatten)]
    input: InputReport,
    results: Vec<WobResult>,
}

#[derive(Debug, Serialize)]
struct Metrics {
    tp_query: u64,
    tp_ref: u64,
    fp: u64,
    #[serde(rename = "fn")]
    fn_: u64,
    precision: Option<f64>,
    recall: Option<f64>,
    f1: Option<f64>,
}

/// All records at the top level, plus multi-exon and mono-exon records separately.
#[derive(Debug, Serialize)]
struct WobResult {
    w: usize,
    #[serde(flatten)]
    all: Metrics,
    multi_exon: Metrics,
    mono_exon: Metrics,
}

fn metrics(tp: u64, query_total: u64, ref_total: u64) -> AnyResult<Metrics> {
    if tp > query_total || tp > ref_total {
        bail!(
            "Matched pairs {tp} exceed transcript totals (query {query_total}, reference {ref_total})"
        );
    }
    let precision = ratio(tp, query_total);
    let recall = ratio(tp, ref_total);
    Ok(Metrics {
        tp_query: tp,
        tp_ref: tp,
        fp: query_total - tp,
        fn_: ref_total - tp,
        precision: round6(precision),
        recall: round6(recall),
        f1: round6(f1_score(precision, recall)),
    })
}

/// Reference summary plus one entry per query sample with results for each w.
/// `inputs[0]` is the reference; denominators come from the evaluated GTFs' index summaries.
fn build_report(
    preset: BenchPreset,
    inputs: Vec<BenchInput>,
    matches: &[SampleMatches],
) -> AnyResult<BenchReport> {
    let totals = |input: &BenchInput| {
        let stats = input.eval_stats();
        (
            stats.transcript_count,
            stats.multi_exon_tx_cnt,
            stats.mono_exon_tx_cnt,
        )
    };
    let report = |sample_idx: usize, input: BenchInput| InputReport {
        total: input.eval_stats().transcript_count,
        sample_id: sample_label(sample_idx),
        input_path: input.raw_path,
        gtf_stats: input.raw_stats,
        collapsed: input.collapsed,
    };
    let mut inputs = inputs.into_iter();
    let reference = inputs.next().context("Missing reference input")?;
    let (ref_all, ref_multi, ref_mono) = totals(&reference);
    let mut queries = Vec::new();
    for (sample_idx, input) in (1..).zip(inputs) {
        let (q_all, q_multi, q_mono) = totals(&input);
        let sample_matches = &matches[sample_idx];
        let mono = sample_matches.mono;
        let mut results = Vec::new();
        for (w, &multi) in sample_matches.multi.iter().enumerate() {
            let context = || format!("{} at w={w}", input.raw_path.display());
            results.push(WobResult {
                w,
                all: metrics(multi + mono, q_all, ref_all).with_context(context)?,
                multi_exon: metrics(multi, q_multi, ref_multi).with_context(context)?,
                mono_exon: metrics(mono, q_mono, ref_mono).with_context(context)?,
            });
        }
        queries.push(QueryReport {
            input: report(sample_idx, input),
            results,
        });
    }
    Ok(BenchReport {
        preset,
        reference: report(0, reference),
        query_count: queries.len(),
        queries,
    })
}

/// results_summary.json: results_all.json without the index summaries, plus the evaluated GTF
/// and the totals behind each ratio.
#[derive(Debug, Serialize)]
struct SummaryReport<'a> {
    preset: BenchPreset,
    /// "input" (default) or "collapsed" (junction): which GTFs the metrics are computed on.
    evaluated_on: &'static str,
    params: SummaryParams,
    reference: SummaryInput<'a>,
    query_count: usize,
    queries: Vec<SummaryQuery<'a>>,
}

/// Thresholds that change the metrics; TSS/TES are null for junction. Full values: params.json.
#[derive(Debug, Serialize)]
struct SummaryParams {
    wob: u32,
    tss_wob: Option<u32>,
    tes_wob: Option<u32>,
    mono_ovlp: f64,
}

/// `input_path`/`input_total`: command-line GTF; `evaluated_path`/`total` and the class totals:
/// GTF the metrics are computed on (the same file for default).
#[derive(Debug, Serialize)]
struct SummaryInput<'a> {
    sample_id: &'a str,
    input_path: &'a Path,
    evaluated_path: &'a Path,
    input_total: u64,
    total: u64,
    multi_exon_total: u64,
    mono_exon_total: u64,
}

#[derive(Debug, Serialize)]
struct SummaryQuery<'a> {
    #[serde(flatten)]
    input: SummaryInput<'a>,
    results: &'a [WobResult],
}

impl InputReport {
    fn summary(&self) -> SummaryInput<'_> {
        let (evaluated_path, eval_stats) = match &self.collapsed {
            Some(c) => (c.path.as_path(), &c.gtf_stats),
            None => (self.input_path.as_path(), &self.gtf_stats),
        };
        SummaryInput {
            sample_id: &self.sample_id,
            input_path: &self.input_path,
            evaluated_path,
            input_total: self.gtf_stats.transcript_count,
            total: self.total,
            multi_exon_total: eval_stats.multi_exon_tx_cnt,
            mono_exon_total: eval_stats.mono_exon_tx_cnt,
        }
    }
}

fn build_summary<'a>(report: &'a BenchReport, args: &BenchArgs) -> SummaryReport<'a> {
    let end_wobs = args.end_wobs();
    SummaryReport {
        preset: report.preset,
        evaluated_on: match report.preset {
            BenchPreset::Default => "input",
            BenchPreset::Junction => "collapsed",
        },
        params: SummaryParams {
            wob: args.wob,
            tss_wob: end_wobs.map(|(tss, _)| tss),
            tes_wob: end_wobs.map(|(_, tes)| tes),
            mono_ovlp: args.mono_ovlp,
        },
        reference: report.reference.summary(),
        query_count: report.query_count,
        queries: report
            .queries
            .iter()
            .map(|q| SummaryQuery {
                input: q.input.summary(),
                results: &q.results,
            })
            .collect(),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TxDiffs {
    tss_diff: u32,
    tes_diff: u32,
    jc_diff_max: u32,
}

/// pure strcture of a transcript
/// used for comparison in a quick way with args
pub struct TxStructure {
    pub sample_idx: usize,
    pub orig_txid: String,
    pub orig_geneid: String,
    pub tss: u32,
    pub tes: u32,
    pub junctions: Vec<u32>,
}

impl TxStructure {
    pub fn is_ref_tx(&self) -> bool {
        self.sample_idx == 0
    }

    /// Exons in genomic order, rebuilt from the terminal coordinates and junctions.
    pub fn exons(&self) -> Vec<(u32, u32)> {
        let bounds: Vec<u32> = std::iter::once(self.tss.min(self.tes))
            .chain(self.junctions.iter().copied())
            .chain(std::iter::once(self.tss.max(self.tes)))
            .collect();
        bounds
            .chunks_exact(2)
            .map(|pair| (pair[0], pair[1]))
            .collect()
    }

    /// Both records come from the same merged group, so junction vectors have equal length.
    pub fn compare(&self, other: &Self) -> TxDiffs {
        TxDiffs {
            tss_diff: self.tss.abs_diff(other.tss),
            tes_diff: self.tes.abs_diff(other.tes),
            jc_diff_max: self
                .junctions
                .iter()
                .zip(&other.junctions)
                .map(|(x, y)| x.abs_diff(*y))
                .max()
                .unwrap_or(0),
        }
    }
}

/// Signed (TSS, TES) differences of `query` minus `reference` along the transcription direction.
/// Positive TSS: query starts downstream (shorter 5'). Positive TES: query ends downstream
/// (longer 3'). Unknown strand follows genomic order, like the TSS/TES convention in TxReader.
/// Matching still uses absolute diffs; this is only for reporting.
// ponytail: not called yet; reserved for the per-transcript output.
#[cfg_attr(not(test), allow(dead_code))]
fn signed_terminal_diffs(
    query: &TxStructure,
    reference: &TxStructure,
    strand: Strand,
) -> (i64, i64) {
    let diff = |q: u32, r: u32| i64::from(q) - i64::from(r);
    let (tss, tes) = (
        diff(query.tss, reference.tss),
        diff(query.tes, reference.tes),
    );
    if strand == Strand::Minus {
        (-tss, -tes)
    } else {
        (tss, tes)
    }
}

pub struct TxGroup {
    chr: String,    // chr from the first column
    tx_id: String,  // merged transcript ID
    strand: Strand, //  strand taken from the record
    tx_structures: Vec<TxStructure>,
    ref_tx_cnt: usize,
    n_exons: usize,
}

/// read merged gtf output from merge function and build TxGroup for each merged transcirpt.
pub struct TxReader {
    reader: BufReader<GzDecoder<BufReader<File>>>,
    line_no: usize,
    done: bool,
    headerview: String,                // header string directly take from the GTF
    sample_list: Vec<(usize, String)>, // One-based sample ID and input file name.
}

impl TxReader {
    pub fn new(path: &Path) -> Result<Self, BenchError> {
        let file = File::open(path).map_err(|source| BenchError::OpenGTF {
            path: path.to_path_buf(),
            source,
        })?;
        let mut tx_reader = Self {
            reader: BufReader::new(GzDecoder::new(BufReader::new(file))),
            line_no: 0,
            done: false,
            headerview: String::new(),
            sample_list: Vec::new(),
        };
        let mut line = String::new();
        // Peek so the first transcript record remains available to read_group.
        while matches!(
            tx_reader.reader.fill_buf()?.first(),
            Some(b'#' | b'\n' | b'\r')
        ) {
            line.clear();
            tx_reader.reader.read_line(&mut line)?;
            tx_reader.line_no += 1;
            if !line.starts_with('#') {
                continue;
            }
            tx_reader.headerview.push_str(&line);
            if line.starts_with("##ISOM <SAMPLE>") {
                let id =
                    parse_header_value(&line, "id").ok_or_else(|| BenchError::MalformedGTF {
                        reason: format!("Missing sample id at line {}", tx_reader.line_no),
                    })?;
                let sample_id = id
                    .strip_prefix('S')
                    .and_then(|value| value.parse::<usize>().ok())
                    .filter(|&value| value > 0)
                    .ok_or_else(|| BenchError::MalformedGTF {
                        reason: format!("Invalid sample id {id} at line {}", tx_reader.line_no),
                    })?;
                let input =
                    parse_header_value(&line, "input").ok_or_else(|| BenchError::MalformedGTF {
                        reason: format!("Missing sample input at line {}", tx_reader.line_no),
                    })?;
                let file_name = Path::new(&input)
                    .file_name()
                    .map(|name| name.to_string_lossy().to_string())
                    .unwrap_or(input);
                tx_reader.sample_list.push((sample_id, file_name));
            }
        }
        Ok(tx_reader)
    }

    fn read_group(&mut self) -> Result<Option<TxGroup>, BenchError> {
        let mut line = String::new();
        loop {
            line.clear();
            if self.reader.read_line(&mut line)? == 0 {
                return Ok(None);
            }
            self.line_no += 1;
            if !line.trim().is_empty() && !line.starts_with('#') {
                break;
            }
        }

        let record = parse_record(line.trim_end())?;
        if record.feature != "transcript" {
            return Err(BenchError::MalformedGTF {
                reason: "Expected transcript record".to_string(),
            });
        }
        let tx_id = attribute(record.attributes, "transcript_id")
            .ok_or_else(|| BenchError::MalformedGTF {
                reason: "Missing transcript_id".to_string(),
            })?
            .to_string();
        let sources = attribute(record.attributes, "ISOM_SRC")
            .ok_or_else(|| BenchError::MalformedGTF {
                reason: "Missing ISOM_SRC; merged GTF must be generated without --chop".to_string(),
            })?
            .to_string();
        let n_exons = attribute(record.attributes, "ISOM_EXONS")
            .ok_or_else(|| BenchError::MalformedGTF {
                reason: "Missing ISOM_EXONS".to_string(),
            })?
            .parse::<usize>()?;
        if n_exons == 0 {
            return Err(BenchError::MalformedGTF {
                reason: format!("Transcript {tx_id} has no exon"),
            });
        }
        let chrom = record.seqid.to_string();
        let strand = record.strand;

        // ponytail: merge output has contiguous exon blocks; arbitrary GTF needs grouping first.
        let mut exons = Vec::new();
        while exons.len() < n_exons {
            line.clear();
            if self.reader.read_line(&mut line)? == 0 {
                return Err(BenchError::MalformedGTF {
                    reason: format!("Incomplete exon block for {tx_id}"),
                });
            }
            self.line_no += 1;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let exon = parse_record(line.trim_end())?;
            if exon.feature != "exon"
                || exon.seqid != chrom
                || exon.strand != strand
                || attribute(exon.attributes, "transcript_id") != Some(tx_id.as_str())
            {
                return Err(BenchError::MalformedGTF {
                    reason: format!("Unexpected record in exon block for {tx_id}"),
                });
            }
            if exon.start == 0 || exon.start > exon.end {
                return Err(BenchError::MalformedGTF {
                    reason: format!("Invalid exon coordinates for {tx_id}"),
                });
            }
            if let Some(&(_, end)) = exons.last()
                && end >= exon.start
            {
                return Err(BenchError::MalformedGTF {
                    reason: format!("Unsorted or overlapping exons for {tx_id}"),
                });
            }
            exons.push((exon.start, exon.end));
        }

        let mut tx_structures = Vec::new();
        let mut ref_tx_cnt = 0;
        for source in sources.split('|') {
            // S#:tx_id:gene_id:start:end:tx_type:donor_diff:acceptor_diff:junction_diffs
            // IDs are not escaped, so split the six fixed fields from the right.
            let malformed = || BenchError::MalformedGTF {
                reason: format!("Malformed ISOM_SRC record: {source}"),
            };
            let (sample, rest) = source.split_once(':').ok_or_else(malformed)?;
            let parts: Vec<&str> = rest.rsplitn(7, ':').collect();
            if parts.len() != 7 {
                return Err(malformed());
            }
            // A gene_id containing ':' is ambiguous: the split keeps transcript IDs intact.
            let (orig_txid, orig_geneid) = parts[6].rsplit_once(':').ok_or_else(malformed)?;
            let file_id = sample
                .strip_prefix('S')
                .ok_or_else(|| BenchError::MalformedGTF {
                    reason: "Missing source sample prefix S".to_string(),
                })?
                .parse::<usize>()?;
            let sample_idx = file_id
                .checked_sub(1)
                .ok_or_else(|| BenchError::MalformedGTF {
                    reason: "Source sample index must start at S1".to_string(),
                })?;
            let start = parts[5].parse::<u32>()?;
            let end = parts[4].parse::<u32>()?;
            let mut source_exons = exons.clone();
            source_exons[0].0 = start;
            source_exons[n_exons - 1].1 = end;
            for (junction_no, left_diff, right_diff) in parse_junction_diffs(parts[0])? {
                if junction_no == 0 || junction_no >= n_exons {
                    return Err(BenchError::MalformedGTF {
                        reason: format!("Invalid junction index in {source}"),
                    });
                }
                source_exons[junction_no - 1].1 =
                    shift_coordinate(exons[junction_no - 1].1, left_diff)?;
                source_exons[junction_no].0 = shift_coordinate(exons[junction_no].0, right_diff)?;
            }
            if !(source_exons
                .iter()
                .all(|&(left, right)| left > 0 && left <= right)
                && source_exons.windows(2).all(|pair| pair[0].1 < pair[1].0))
            {
                return Err(BenchError::MalformedGTF {
                    reason: format!("Invalid reconstructed exons in {source}"),
                });
            }
            let (tss, tes) = if strand == Strand::Minus {
                (end, start)
            } else {
                (start, end)
            };
            ref_tx_cnt += usize::from(sample_idx == 0);
            tx_structures.push(TxStructure {
                sample_idx,
                orig_txid: orig_txid.to_string(),
                orig_geneid: orig_geneid.to_string(),
                tss,
                tes,
                // Genomic order: exon end, next exon start, repeated for each junction.
                junctions: source_exons
                    .windows(2)
                    .flat_map(|pair| [pair[0].1, pair[1].0])
                    .collect(),
            });
        }
        Ok(Some(TxGroup {
            chr: chrom,
            tx_id,
            strand,
            tx_structures,
            ref_tx_cnt,
            n_exons,
        }))
    }
}

impl Iterator for TxReader {
    type Item = Result<TxGroup, BenchError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        match self.read_group() {
            Ok(Some(group)) => Some(Ok(group)),
            Ok(None) => {
                self.done = true;
                None
            }
            Err(err) => {
                self.done = true;
                Some(Err(BenchError::ReadGTF {
                    line_no: self.line_no,
                    source: Box::new(err),
                }))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Cli, Commands};
    use clap::Parser;
    use flate2::{Compression, write::GzEncoder};
    use std::io::Write;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn bench_args(extra: &[&str]) -> BenchArgs {
        let mut cli = vec![
            "isomatch", "bench", "-g", "ref.gtf", "-o", "out", "-r", "ref.fa",
        ];
        cli.extend_from_slice(extra);
        cli.push("query.gtf");
        let Commands::Bench(args) = Cli::try_parse_from(cli).unwrap().command else {
            panic!("expected bench arguments");
        };
        args
    }

    /// Build a source record from exons the same way TxReader does.
    fn tx(sample_idx: usize, id: &str, strand: Strand, exons: &[(u32, u32)]) -> TxStructure {
        let (start, end) = (exons[0].0, exons[exons.len() - 1].1);
        let (tss, tes) = if strand == Strand::Minus {
            (end, start)
        } else {
            (start, end)
        };
        TxStructure {
            sample_idx,
            orig_txid: id.to_string(),
            orig_geneid: format!("g_{id}"),
            tss,
            tes,
            junctions: exons
                .windows(2)
                .flat_map(|pair| [pair[0].1, pair[1].0])
                .collect(),
        }
    }

    fn group(tx_structures: Vec<TxStructure>) -> TxGroup {
        let n_exons = tx_structures[0].junctions.len() / 2 + 1;
        TxGroup {
            chr: "chr1".to_string(),
            tx_id: "ISOMT_1".to_string(),
            strand: Strand::Plus,
            ref_tx_cnt: tx_structures.iter().filter(|tx| tx.is_ref_tx()).count(),
            tx_structures,
            n_exons,
        }
    }

    /// Pairs per sample for each w (multi-exon pairs at w plus mono-exon pairs).
    fn run_group(args: &BenchArgs, txs: Vec<TxStructure>, n_samples: usize) -> Vec<Vec<u64>> {
        let mut matches = vec![
            SampleMatches {
                multi: vec![0; args.wob as usize + 1],
                mono: 0,
            };
            n_samples
        ];
        match_group(&group(txs), args, &mut matches);
        matches
            .iter()
            .map(|m| m.multi.iter().map(|&multi| multi + m.mono).collect())
            .collect()
    }

    #[test]
    fn signed_terminal_diffs_follow_transcription_direction() {
        let exons_ref: &[(u32, u32)] = &[(100, 200), (300, 400)];
        // Query starts 5 bp right of and ends 10 bp right of the reference in genomic order.
        let exons_q: &[(u32, u32)] = &[(105, 200), (300, 410)];
        for (strand, expected) in [
            // +: 5' end shorter (TSS downstream), 3' end longer (TES downstream).
            (Strand::Plus, (5, 10)),
            // -: TSS is the right end and moved 10 bp upstream (longer 5'); TES is the left end
            // and moved 5 bp upstream (shorter 3').
            (Strand::Minus, (-10, -5)),
            (Strand::Unknown, (5, 10)),
        ] {
            let r = tx(0, "r", strand, exons_ref);
            let q = tx(1, "q", strand, exons_q);
            assert_eq!(
                signed_terminal_diffs(&q, &r, strand),
                expected,
                "{strand:?}"
            );
            assert_eq!(
                signed_terminal_diffs(&r, &q, strand),
                (-expected.0, -expected.1)
            );
        }
    }

    #[test]
    fn candidate_thresholds_are_inclusive() {
        let args = bench_args(&["-w", "2", "-s", "5", "-e", "5"]);
        let plus = Strand::Plus;
        let r = tx(0, "r", plus, &[(100, 200), (300, 400), (500, 600)]);
        let ok = tx(1, "q", plus, &[(105, 202), (298, 400), (500, 595)]);
        assert_eq!(
            candidate_diffs(&ok, &r, &args),
            Some(TxDiffs {
                tss_diff: 5,
                tes_diff: 5,
                jc_diff_max: 2
            })
        );
        for bad in [
            tx(1, "q", plus, &[(94, 200), (300, 400), (500, 600)]),
            tx(1, "q", plus, &[(100, 200), (300, 400), (500, 606)]),
            tx(1, "q", plus, &[(100, 203), (300, 400), (500, 600)]),
        ] {
            assert_eq!(candidate_diffs(&bad, &r, &args), None);
        }
        // Minus strand: TSS is the right end; diffs use abs_diff in both directions.
        let minus = Strand::Minus;
        let r = tx(0, "r", minus, &[(100, 200), (300, 400)]);
        let q = tx(1, "q", minus, &[(110, 200), (300, 395)]);
        assert_eq!(
            candidate_diffs(&q, &r, &args),
            None,
            "TES (left end) differs by 10"
        );
        let q = tx(1, "q", minus, &[(96, 200), (300, 405)]);
        assert_eq!(
            candidate_diffs(&q, &r, &args),
            Some(TxDiffs {
                tss_diff: 5,
                tes_diff: 4,
                jc_diff_max: 0
            })
        );
    }

    #[test]
    fn mono_candidates_need_overlap_and_reciprocal_fraction() {
        let args = bench_args(&["-s", "1000", "-e", "1000", "--mono-ovlp", "0.9"]);
        let plus = Strand::Plus;
        let r = tx(0, "r", plus, &[(1, 100)]);
        assert!(candidate_diffs(&tx(1, "q", plus, &[(1, 90)]), &r, &args).is_some());
        assert!(candidate_diffs(&tx(1, "q", plus, &[(1, 89)]), &r, &args).is_none());
        assert!(candidate_diffs(&tx(1, "q", plus, &[(1, 111)]), &r, &args).is_some());
        assert!(candidate_diffs(&tx(1, "q", plus, &[(1, 112)]), &r, &args).is_none());
        let one_bp = tx(0, "r", plus, &[(50, 50)]);
        assert!(candidate_diffs(&tx(1, "q", plus, &[(50, 50)]), &one_bp, &args).is_some());
        assert!(candidate_diffs(&tx(1, "q", plus, &[(51, 51)]), &one_bp, &args).is_none());
        // Same merged group via transitive links but no direct overlap: rejected; a shared end counts.
        let loose = bench_args(&["-s", "1000", "-e", "1000", "--mono-ovlp", "0.005"]);
        assert!(candidate_diffs(&tx(1, "q", plus, &[(101, 200)]), &r, &loose).is_none());
        assert!(candidate_diffs(&tx(1, "q", plus, &[(100, 200)]), &r, &loose).is_some());
    }

    #[test]
    fn greedy_matching_follows_priority_per_threshold() {
        let args = bench_args(&["-w", "1"]);
        let p = Strand::Plus;
        // Priority Q1-R1 (jc 0), Q1-R2 (jc 1), Q2-R1 (jc 1, larger TSS diff): one pair.
        let txs = vec![
            tx(0, "R1", p, &[(100, 200), (300, 400)]),
            tx(0, "R2", p, &[(100, 201), (300, 400)]),
            tx(1, "Q1", p, &[(100, 200), (300, 400)]),
            tx(1, "Q2", p, &[(110, 199), (300, 400)]),
        ];
        assert_eq!(run_group(&args, txs, 2)[1], [1, 1]);
        // Q2 loses R1 but takes the alternative R2 once w allows it.
        let txs = vec![
            tx(0, "R1", p, &[(100, 200), (300, 400)]),
            tx(0, "R2", p, &[(100, 201), (300, 400)]),
            tx(1, "Q1", p, &[(100, 200), (300, 400)]),
            tx(1, "Q2", p, &[(110, 200), (300, 400)]),
        ];
        assert_eq!(run_group(&args, txs, 2)[1], [1, 2]);
    }

    #[test]
    fn ties_break_by_transcript_id() {
        let args = bench_args(&["-s", "10", "-e", "10"]);
        let p = Strand::Plus;
        // Both queries differ from R1 by TSS 5, a tie. Only `alt` also fits R2 (TSS 5, TES 5).
        let txs = |alt: &str, other: &str| {
            vec![
                tx(0, "R1", p, &[(100, 200), (300, 400)]),
                tx(0, "R2", p, &[(90, 200), (300, 405)]),
                tx(1, alt, p, &[(95, 200), (300, 400)]),
                tx(1, other, p, &[(105, 200), (300, 400)]),
            ]
        };
        // "a" wins the tie on R1, so its alternative is wasted and "b" stays unmatched.
        assert_eq!(run_group(&args, txs("a", "b"), 2)[1], [1]);
        // With IDs swapped, the query without an alternative takes R1 and "b" takes R2.
        assert_eq!(run_group(&args, txs("b", "a"), 2)[1], [2]);
    }

    #[test]
    fn samples_are_matched_independently_and_mono_ignores_w() {
        let args = bench_args(&["-w", "2"]);
        let p = Strand::Plus;
        let txs = vec![
            tx(0, "R", p, &[(100, 200)]),
            tx(1, "Q", p, &[(101, 200)]),
            tx(2, "Q", p, &[(100, 199)]),
            tx(2, "Q_dup", p, &[(100, 199)]),
        ];
        let matches = run_group(&args, txs, 3);
        assert_eq!(matches[1], [1, 1, 1]);
        assert_eq!(matches[2], [1, 1, 1]);
        // A group without reference records adds nothing.
        let matches = run_group(&args, vec![tx(1, "Q", p, &[(1, 10)])], 2);
        assert_eq!(matches[1], [0, 0, 0]);
    }

    #[test]
    fn junction_preset_ignores_ends_except_mono_overlap() {
        let args = bench_args(&["--preset", "junction", "-w", "1", "--mono-ovlp", "0.9"]);
        assert_eq!(args.end_wobs(), None);
        assert_eq!(bench_args(&[]).end_wobs(), Some((50, 50)));
        assert_eq!(bench_args(&["-s", "7"]).end_wobs(), Some((7, 50)));
        let p = Strand::Plus;
        let r = tx(0, "r", p, &[(100, 200), (300, 400), (500, 600)]);
        // Ends far beyond the default thresholds; diffs stay the real values.
        let q = tx(1, "q", p, &[(10, 201), (300, 400), (500, 900)]);
        assert_eq!(
            candidate_diffs(&q, &r, &args),
            Some(TxDiffs {
                tss_diff: 90,
                tes_diff: 300,
                jc_diff_max: 1
            })
        );
        assert_eq!(candidate_diffs(&q, &r, &bench_args(&["-w", "1"])), None);
        let q = tx(1, "q", p, &[(100, 202), (300, 400), (500, 600)]);
        assert_eq!(candidate_diffs(&q, &r, &args), None, "junction diff 2 > w");
        // Mono: only reciprocal overlap, no TSS/TES threshold.
        let r = tx(0, "r", p, &[(1, 1000)]);
        assert!(candidate_diffs(&tx(1, "q", p, &[(81, 1000)]), &r, &args).is_some());
        assert!(candidate_diffs(&tx(1, "q", p, &[(102, 1000)]), &r, &args).is_none());
    }

    #[test]
    fn junction_preset_breaks_multi_exon_ties_by_id_not_ends() {
        let p = Strand::Plus;
        // Both queries have the reference chain; "b" has closer ends, "a" the smaller ID.
        let txs = || {
            vec![
                tx(0, "R", p, &[(100, 200), (300, 400)]),
                tx(1, "a", p, &[(150, 200), (300, 450)]),
                tx(1, "b", p, &[(100, 200), (300, 400)]),
            ]
        };
        let pair = |args: &BenchArgs| {
            let group = group(txs());
            let mut matches = vec![SampleMatches::default(); 2];
            matches[1].multi = vec![0];
            let pairs = match_group(&group, args, &mut matches);
            group.tx_structures[pairs[0].0].orig_txid.clone()
        };
        assert_eq!(pair(&bench_args(&["--preset", "junction"])), "a");
        assert_eq!(pair(&bench_args(&["-s", "100", "-e", "100"])), "b");
        // Mono-exon groups keep the end-based order under junction.
        let mono = vec![
            tx(0, "R", p, &[(100, 400)]),
            tx(1, "a", p, &[(110, 400)]),
            tx(1, "b", p, &[(100, 400)]),
        ];
        let group = group(mono);
        let mut matches = vec![SampleMatches::default(); 2];
        let pairs = match_group(&group, &bench_args(&["--preset", "junction"]), &mut matches);
        assert_eq!(group.tx_structures[pairs[0].0].orig_txid, "b");
        assert_eq!(matches[1].mono, 1);
    }

    #[test]
    fn prepare_outdir_only_clears_marked_bench_output() -> AnyResult<()> {
        let dir = temp_dir("outdir")?;
        let marker = |out: &Path| out.join(BENCH_MARKER).is_file();
        let outside = dir.join("input.gtf");
        std::fs::write(&outside, "keep")?;

        // Missing directory: created with the marker, with or without --overwrite.
        let out = dir.join("new");
        prepare_outdir(&out, false, &[], &[])?;
        assert!(marker(&out));
        // Empty directory: reused.
        let out = dir.join("empty");
        std::fs::create_dir_all(&out)?;
        prepare_outdir(&out, false, &[], &[])?;
        assert!(marker(&out));

        // Non-empty without marker: refused even with --overwrite; nothing deleted.
        let out = dir.join("foreign");
        std::fs::create_dir_all(&out)?;
        std::fs::write(out.join("notes.txt"), "user file")?;
        for overwrite in [false, true] {
            assert!(prepare_outdir(&out, overwrite, &[], &[]).is_err());
        }
        assert!(out.join("notes.txt").exists() && !marker(&out));

        // Marked bench output: refused without --overwrite, cleared with it.
        let out = dir.join("bench");
        prepare_outdir(&out, false, &[], &[])?;
        std::fs::create_dir_all(out.join("Q1"))?;
        std::fs::write(out.join("Q1").join("fp.gtf.gz"), "old")?;
        #[cfg(unix)]
        std::os::unix::fs::symlink(&outside, out.join("link"))?;
        assert!(prepare_outdir(&out, false, &[], &[]).is_err());
        // Protected paths inside it (an input, the working directory) or a forbidden target block it.
        let inside = out.join("query.gtf");
        std::fs::write(&inside, "query")?;
        let err = prepare_outdir(&out, true, &[outside.clone(), inside.clone()], &[]).unwrap_err();
        assert!(err.to_string().contains("query.gtf"));
        assert!(prepare_outdir(&out, true, &[out.join("Q1")], &[]).is_err());
        std::fs::remove_file(&inside)?;
        assert!(prepare_outdir(&out, true, &[], std::slice::from_ref(&out)).is_err());
        assert!(out.join("Q1").join("fp.gtf.gz").exists());

        prepare_outdir(
            &out,
            true,
            &[outside.clone(), dir.join("missing.gtf")],
            std::slice::from_ref(&dir),
        )?;
        let left: Vec<_> = std::fs::read_dir(&out)?
            .map(|e| e.unwrap().file_name())
            .collect();
        assert_eq!(left, [BENCH_MARKER]);
        assert_eq!(
            std::fs::read_to_string(&outside)?,
            "keep",
            "symlink target untouched"
        );
        std::fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn metrics_and_report() -> AnyResult<()> {
        assert_eq!(ratio(1, 0), None);
        assert_eq!(f1_score(None, Some(1.0)), None);
        assert_eq!(f1_score(Some(0.0), Some(0.0)), Some(0.0));
        assert_eq!(round6(f1_score(Some(0.5), Some(1.0))), Some(0.666667));
        assert_eq!(round6(None), None);
        assert!(metrics(3, 2, 5).is_err());
        assert!(metrics(3, 5, 2).is_err());

        // (multi, mono) transcript counts.
        let stats = |(multi, mono): (u64, u64)| IndexStats {
            transcript_count: multi + mono,
            multi_exon_tx_cnt: multi,
            mono_exon_tx_cnt: mono,
            ..Default::default()
        };
        let input = |path: &str, counts, collapsed: Option<(u64, u64)>| BenchInput {
            raw_path: PathBuf::from(path),
            raw_stats: stats(counts),
            collapsed: collapsed.map(|c| CollapsedInput {
                path: PathBuf::from(format!("c_{path}")),
                gtf_stats: stats(c),
            }),
        };
        let inputs = vec![
            input("ref.gtf", (1, 1), None),
            input("q1.gtf", (2, 2), None),
            input("q2.gtf", (2, 0), None),
        ];
        let matches = [
            SampleMatches::default(),
            SampleMatches {
                multi: vec![0, 1],
                mono: 1,
            },
            SampleMatches {
                multi: vec![0, 0],
                mono: 0,
            },
        ];
        let json = serde_json::to_value(build_report(BenchPreset::Default, inputs, &matches)?)?;
        assert_eq!(json["preset"], "default");
        assert_eq!(json["query_count"], 2);
        assert_eq!(json["reference"]["sample_id"], "REF");
        assert_eq!(json["reference"]["input_path"], "ref.gtf");
        assert_eq!(json["reference"]["total"], 2);
        assert_eq!(json["reference"]["gtf_stats"]["mono_exon_tx_cnt"], 1);
        assert!(json["reference"]["collapsed"].is_null());
        let q1 = &json["queries"][0];
        assert_eq!(
            (&q1["sample_id"], &q1["input_path"], &q1["total"]),
            (&"Q1".into(), &"q1.gtf".into(), &4.into())
        );
        assert_eq!(
            q1["results"][1],
            serde_json::json!({
                "w": 1, "tp_query": 2, "tp_ref": 2, "fp": 2, "fn": 0,
                "precision": 0.5, "recall": 1.0, "f1": 0.666667,
                "multi_exon": {"tp_query": 1, "tp_ref": 1, "fp": 1, "fn": 0,
                               "precision": 0.5, "recall": 1.0, "f1": 0.666667},
                "mono_exon": {"tp_query": 1, "tp_ref": 1, "fp": 1, "fn": 0,
                              "precision": 0.5, "recall": 1.0, "f1": 0.666667}
            })
        );
        assert_eq!(q1["results"][0]["tp_query"], 1);
        assert_eq!(q1["results"][0]["multi_exon"]["tp_query"], 0);
        let q2 = &json["queries"][1]["results"][0];
        assert_eq!(q2["f1"], 0.0);
        assert!(q2["mono_exon"]["precision"].is_null());

        // Junction: totals come from the collapsed GTFs; gtf_stats stays the raw one.
        let inputs = vec![
            input("ref.gtf", (5, 1), Some((3, 1))),
            input("q1.gtf", (4, 0), Some((2, 0))),
        ];
        let matches = [
            SampleMatches::default(),
            SampleMatches {
                multi: vec![2],
                mono: 0,
            },
        ];
        let json = serde_json::to_value(build_report(BenchPreset::Junction, inputs, &matches)?)?;
        assert_eq!(json["preset"], "junction");
        assert_eq!(json["reference"]["total"], 4);
        assert_eq!(json["reference"]["gtf_stats"]["transcript_count"], 6);
        assert_eq!(json["reference"]["collapsed"]["path"], "c_ref.gtf");
        assert_eq!(
            json["reference"]["collapsed"]["gtf_stats"]["transcript_count"],
            4
        );
        assert_eq!(json["queries"][0]["results"][0]["precision"], 1.0);
        assert_eq!(json["queries"][0]["results"][0]["recall"], 0.5);

        let inputs = vec![
            input("ref.gtf", (1, 0), None),
            input("q1.gtf", (2, 0), None),
        ];
        let too_many = [
            SampleMatches::default(),
            SampleMatches {
                multi: vec![2],
                mono: 0,
            },
        ];
        assert!(build_report(BenchPreset::Default, inputs, &too_many).is_err());
        Ok(())
    }

    fn temp_dir(tag: &str) -> AnyResult<PathBuf> {
        let dir = std::env::temp_dir().join(format!(
            "isomatch-bench-{tag}-{}-{}",
            std::process::id(),
            SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
        ));
        std::fs::create_dir_all(&dir)?;
        Ok(dir)
    }

    type TxDef = (&'static str, char, &'static [(u32, u32)]);

    fn write_gtf(path: &Path, txs: &[TxDef]) -> AnyResult<()> {
        let mut out = String::new();
        let mut sorted = txs.to_vec();
        sorted.sort_by_key(|(_, _, exons)| (exons[0].0, exons[exons.len() - 1].1));
        for (id, strand, exons) in sorted {
            let (start, end) = (exons[0].0, exons[exons.len() - 1].1);
            let attr = format!("gene_id \"g_{id}\"; transcript_id \"{id}\";");
            out.push_str(&format!(
                "chr1\tt\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t{attr}\n"
            ));
            for (s, e) in exons.iter() {
                out.push_str(&format!(
                    "chr1\tt\texon\t{s}\t{e}\t.\t{strand}\t.\t{attr}\n"
                ));
            }
        }
        std::fs::write(path, out)?;
        Ok(())
    }

    /// Reference result: compare every query with every reference genome-wide, no merge groups.
    fn brute_force(args: &BenchArgs, refs: &[TxDef], queries: &[TxDef]) -> Vec<u64> {
        let to_tx = |sample_idx, (id, strand, exons): &TxDef| {
            let strand = if *strand == '-' {
                Strand::Minus
            } else {
                Strand::Plus
            };
            (strand, tx(sample_idx, id, strand, exons))
        };
        let refs: Vec<_> = refs.iter().map(|t| to_tx(0, t)).collect();
        let qs: Vec<_> = queries.iter().map(|t| to_tx(1, t)).collect();
        let mut cands = Vec::new();
        for (qi, (qs_strand, q)) in qs.iter().enumerate() {
            for (ri, (rs_strand, r)) in refs.iter().enumerate() {
                if qs_strand == rs_strand
                    && q.junctions.len() == r.junctions.len()
                    && let Some(d) = candidate_diffs(q, r, args)
                {
                    cands.push((d, qi, ri));
                }
            }
        }
        cands.sort_by(|(a, aq, ar), (b, bq, br)| {
            (
                a.jc_diff_max,
                a.tss_diff,
                a.tes_diff,
                &qs[*aq].1.orig_txid,
                &refs[*ar].1.orig_txid,
            )
                .cmp(&(
                    b.jc_diff_max,
                    b.tss_diff,
                    b.tes_diff,
                    &qs[*bq].1.orig_txid,
                    &refs[*br].1.orig_txid,
                ))
        });
        (0..=args.wob as usize)
            .map(|w| {
                let n = cands.partition_point(|(d, _, _)| d.jc_diff_max as usize <= w);
                greedy_match(&cands[..n], qs.len(), refs.len()).len() as u64
            })
            .collect()
    }

    /// Flatten the queries in results_all.json into rows: sample_id, input_path, w, query_total, reference_total,
    /// tp_query, tp_ref, fp, fn, precision, recall, f1.
    /// transcript_id -> (strand, exons) from `<sample_dir>/<status>.gtf.gz`.
    type StatusTxs = BTreeMap<String, (char, Vec<(u32, u32)>)>;

    fn read_status_gtf(sample_dir: &Path, status: &str) -> AnyResult<StatusTxs> {
        let file = File::open(sample_dir.join(format!("{status}.gtf.gz")))?;
        let mut text = String::new();
        std::io::Read::read_to_string(&mut GzDecoder::new(BufReader::new(file)), &mut text)?;
        let mut txs: BTreeMap<String, (char, Vec<(u32, u32)>)> = BTreeMap::new();
        for line in text.lines().filter(|line| !line.starts_with('#')) {
            let record = parse_record(line)?;
            let id = attribute(record.attributes, "transcript_id")
                .unwrap()
                .to_string();
            let entry = txs
                .entry(id)
                .or_insert((char::from(record.strand), Vec::new()));
            if record.feature == "exon" {
                entry.1.push((record.start, record.end));
            }
        }
        Ok(txs)
    }

    fn report_rows(outdir: &Path) -> AnyResult<Vec<Vec<String>>> {
        let json: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(outdir.join("results_all.json"))?)?;
        let mut rows = Vec::new();
        assert_eq!(
            json["query_count"],
            json["queries"].as_array().unwrap().len()
        );
        // total counts the evaluated GTF: the collapsed one under junction.
        let eval_count = |entry: &serde_json::Value| {
            let stats = if entry["collapsed"].is_null() {
                &entry["gtf_stats"]
            } else {
                &entry["collapsed"]["gtf_stats"]
            };
            stats["transcript_count"].clone()
        };
        let reference = &json["reference"];
        assert_eq!(eval_count(reference), reference["total"]);
        for sample in json["queries"].as_array().unwrap() {
            for result in sample["results"].as_array().unwrap() {
                let sample_fields = ["sample_id", "input_path"].map(|k| &sample[k]);
                assert_eq!(eval_count(sample), sample["total"]);
                let totals = [&sample["total"], &json["reference"]["total"]];
                let metrics = [
                    "tp_query",
                    "tp_ref",
                    "fp",
                    "fn",
                    "precision",
                    "recall",
                    "f1",
                ]
                .map(|k| &result[k]);
                let row = sample_fields
                    .into_iter()
                    .chain([&result["w"]])
                    .chain(totals)
                    .chain(metrics)
                    .map(|v| v.as_str().map_or_else(|| v.to_string(), str::to_string))
                    .collect();
                rows.push(row);
            }
        }
        Ok(rows)
    }

    #[test]
    fn bench_end_to_end_matches_brute_force() -> AnyResult<()> {
        let dir = temp_dir("e2e")?;
        let seq: String = "ACGTGTAAGCAGGTCCTTAG".repeat(150);
        std::fs::write(dir.join("ref.fa"), format!(">chr1\n{seq}\n"))?;
        std::fs::write(
            dir.join("ref.fa.fai"),
            format!("chr1\t{}\t6\t{}\t{}\n", seq.len(), seq.len(), seq.len() + 1),
        )?;

        const REF: &[TxDef] = &[
            ("r1", '+', &[(100, 200), (300, 400), (500, 600)]),
            ("r2", '+', &[(100, 200), (303, 400), (500, 600)]),
            ("r3", '-', &[(1000, 1100), (1200, 1300)]),
            ("r4", '+', &[(1500, 1700)]),
            ("r5", '+', &[(2000, 2100), (2200, 2300)]),
            ("r6", '+', &[(2600, 2700)]),
        ];
        const Q1: &[TxDef] = &[
            ("a", '+', &[(105, 200), (301, 400), (500, 590)]),
            ("b", '+', &[(100, 200), (302, 400), (500, 600)]),
            ("c", '+', &[(100, 200), (300, 400), (500, 600)]),
            ("d", '-', &[(1010, 1100), (1202, 1300)]),
            ("e", '+', &[(1510, 1700)]),
            ("f", '+', &[(2000, 2100), (2210, 2300)]),
            ("g", '+', &[(2050, 2100), (2200, 2250)]),
            ("h", '+', &[(2900, 2950)]),
        ];
        const Q2: &[TxDef] = &[
            ("x", '+', &[(100, 200), (300, 400), (500, 600)]),
            ("y", '+', &[(1400, 1520)]),
            ("z", '-', &[(1000, 1100), (1203, 1300)]),
            ("w", '+', &[(100, 200), (301, 400), (500, 600)]),
        ];
        for (name, txs) in [("ref.gtf", REF), ("q1.gtf", Q1), ("q2.gtf", Q2)] {
            write_gtf(&dir.join(name), txs)?;
        }

        let run = |out: &str, queries: &[&str]| -> AnyResult<Vec<Vec<String>>> {
            let ref_gtf = dir.join("ref.gtf").to_string_lossy().to_string();
            let ref_fa = dir.join("ref.fa").to_string_lossy().to_string();
            let outdir = dir.join(out).to_string_lossy().to_string();
            let mut cli = vec![
                "isomatch",
                "bench",
                "-g",
                &ref_gtf,
                "-o",
                &outdir,
                "-r",
                &ref_fa,
                "-w",
                "3",
                "-s",
                "20",
                "-e",
                "20",
                "--mono-ovlp",
                "0.8",
            ];
            let query_paths: Vec<String> = queries
                .iter()
                .map(|q| dir.join(q).to_string_lossy().to_string())
                .collect();
            cli.extend(query_paths.iter().map(String::as_str));
            let Commands::Bench(args) = Cli::try_parse_from(cli).unwrap().command else {
                panic!("expected bench arguments");
            };
            run_bench(args)?;
            report_rows(&dir.join(out))
        };

        let args = bench_args(&["-w", "3", "-s", "20", "-e", "20", "--mono-ovlp", "0.8"]);
        let rows = run("both", &["q1.gtf", "q2.gtf"])?;
        assert_eq!(rows.len(), 2 * 4);
        for (sample, queries, offset) in [("Q1", Q1, 0), ("Q2", Q2, 4)] {
            let expected = brute_force(&args, REF, queries);
            for w in 0..4 {
                let row = &rows[offset + w];
                assert_eq!(
                    (row[0].as_str(), row[2].as_str()),
                    (sample, w.to_string().as_str())
                );
                assert_eq!(row[3], queries.len().to_string());
                assert_eq!(row[4], REF.len().to_string());
                assert_eq!(row[5], expected[w].to_string(), "{sample} w={w}");
            }
        }
        // Brute-force expectations for q1 are not trivial: w changes the result.
        assert_ne!(rows[0][5], rows[3][5]);

        // Status GTF header follows the merge style.
        let file = File::open(dir.join("both").join("Q2").join("fn.gtf.gz"))?;
        let mut text = String::new();
        std::io::Read::read_to_string(&mut GzDecoder::new(BufReader::new(file)), &mut text)?;
        let header: Vec<&str> = text.lines().take_while(|l| l.starts_with('#')).collect();
        assert_eq!(header.len(), 6);
        assert!(header[0].starts_with("##ISOM <VERSION>") && header[0].ends_with("schema=\"1\""));
        assert!(header[1].starts_with(
            "##ISOM <STATUS_BENCH> id=\"fn\"; side=\"reference\"; w=\"0\"; preset=\"default\"; description=\"reference transcripts not matched by any Q2 transcript at w=0\";"
        ));
        let ref_gtf = dir.join("ref.gtf").to_string_lossy().to_string();
        let q2_gtf = dir.join("q2.gtf").to_string_lossy().to_string();
        assert_eq!(
            header[2],
            format!(
                "##ISOM <SAMPLE> id=\"REF\"; input=\"{ref_gtf}\"; role=\"reference\"; total=\"{}\";",
                REF.len()
            )
        );
        assert_eq!(
            header[3],
            format!(
                "##ISOM <SAMPLE> id=\"Q2\"; input=\"{q2_gtf}\"; role=\"query\"; total=\"{}\";",
                Q2.len()
            )
        );
        assert!(header[4].starts_with("##ISOM <FORMAT> ID=\"ISOM_TX_ID\";"));
        assert!(header[5].starts_with("##ISOM <COMMAND> cmd="));

        // Status GTFs at w = 0: each record lands in exactly one file with its original exons.
        for (sample, queries) in [("Q1", Q1), ("Q2", Q2)] {
            let tp_count = brute_force(&args, REF, queries)[0];
            let read = |status: &str| read_status_gtf(&dir.join("both").join(sample), status);
            let (tp_query, fp, tp_ref, fn_) =
                (read("tp-query")?, read("fp")?, read("tp-ref")?, read("fn")?);
            assert_eq!(tp_query.len() as u64, tp_count, "{sample}");
            assert_eq!(tp_ref.len() as u64, tp_count, "{sample}");
            for (defs, a, b) in [(queries, &tp_query, &fp), (REF, &tp_ref, &fn_)] {
                assert_eq!(a.len() + b.len(), defs.len(), "{sample}");
                for (id, strand, exons) in defs {
                    let found = a.get(*id).or_else(|| b.get(*id));
                    assert_eq!(found, Some(&(*strand, exons.to_vec())), "{sample} {id}");
                }
            }
        }

        // Reordering or dropping other queries does not change a sample's metrics.
        let alone = run("alone", &["q1.gtf"])?;
        let reordered = run("reordered", &["q2.gtf", "q1.gtf"])?;
        for w in 0..4 {
            assert_eq!(alone[w][5..], rows[w][5..]);
            assert_eq!(reordered[4 + w][5..], rows[w][5..]);
            assert_eq!(reordered[w][5..], rows[4 + w][5..]);
        }
        std::fs::remove_dir_all(dir)?;
        Ok(())
    }

    /// `isomatch bench` on GTFs in `dir`, run through the real CLI parser; returns results_all.json.
    fn run_cli(
        dir: &Path,
        out: &str,
        ref_gtf: &str,
        queries: &[&str],
        extra: &[&str],
    ) -> AnyResult<serde_json::Value> {
        let path = |name: &str| dir.join(name).to_string_lossy().to_string();
        let (ref_gtf, ref_fa, outdir) = (path(ref_gtf), path("ref.fa"), path(out));
        let mut cli = vec![
            "isomatch", "bench", "-g", &ref_gtf, "-o", &outdir, "-r", &ref_fa,
        ];
        cli.extend_from_slice(extra);
        let query_paths: Vec<String> = queries.iter().map(|q| path(q)).collect();
        cli.extend(query_paths.iter().map(String::as_str));
        let Commands::Bench(args) = Cli::try_parse_from(cli).unwrap().command else {
            panic!("expected bench arguments");
        };
        run_bench(args)?;
        Ok(serde_json::from_str(&std::fs::read_to_string(
            dir.join(out).join("results_all.json"),
        )?)?)
    }

    #[test]
    fn junction_preset_end_to_end() -> AnyResult<()> {
        let dir = temp_dir("junction")?;
        let seq: String = "ACGTGTAAGCAGGTCCTTAG".repeat(150);
        std::fs::write(dir.join("ref.fa"), format!(">chr1\n{seq}\n"))?;
        std::fs::write(
            dir.join("ref.fa.fai"),
            format!("chr1\t{}\t6\t{}\t{}\n", seq.len(), seq.len(), seq.len() + 1),
        )?;
        const REF: &[TxDef] = &[
            ("r1", '+', &[(100, 200), (300, 400), (500, 600)]),
            ("r1b", '+', &[(90, 200), (300, 400), (500, 650)]), // r1 chain, other ends
            ("r2", '+', &[(100, 200), (303, 400), (500, 600)]),
            ("r3", '+', &[(1000, 1100), (1200, 1300), (1400, 1500)]),
            ("r4", '+', &[(1700, 1800)]),
            ("r4dup", '+', &[(1700, 1800)]), // exact mono duplicate
            ("r5", '+', &[(1900, 2000)]),
        ];
        const Q1: &[TxDef] = &[
            ("a", '+', &[(150, 200), (300, 400), (500, 560)]), // r1 chain, ends off by 50/40+
            ("a2", '+', &[(120, 200), (300, 400), (500, 590)]), // a chain
            ("b", '+', &[(100, 200), (301, 400), (500, 600)]), // jc 1 from r1, 2 from r2
            ("c", '+', &[(300, 400), (500, 600)]),             // partial chain of r1
            ("d", '+', &[(1705, 1800)]),
            ("d2", '+', &[(1706, 1800)]), // 1 bp from d: not collapsed
            ("e", '+', &[(1000, 1100), (1200, 1300)]), // partial chain of r3
        ];
        write_gtf(&dir.join("ref.gtf"), REF)?;
        write_gtf(&dir.join("q1.gtf"), Q1)?;

        let json = run_cli(
            &dir,
            "out",
            "ref.gtf",
            &["q1.gtf"],
            &["--preset", "junction", "-w", "3"],
        )?;
        assert_eq!(json["preset"], "junction");
        // REF: r1/r1b and r4/r4dup collapse; Q1: a/a2 collapse, d2 stays.
        let reference = &json["reference"];
        assert_eq!(reference["gtf_stats"]["transcript_count"], 7);
        assert_eq!(reference["total"], 5);
        assert_eq!(reference["collapsed"]["gtf_stats"]["multi_exon_tx_cnt"], 3);
        assert_eq!(reference["collapsed"]["gtf_stats"]["mono_exon_tx_cnt"], 2);
        let q1 = &json["queries"][0];
        assert_eq!(q1["total"], 6);
        let results = q1["results"].as_array().unwrap();
        let tp = |w: usize, class: &str| {
            let r = if class.is_empty() {
                &results[w]
            } else {
                &results[w][class]
            };
            r["tp_query"].as_u64().unwrap()
        };
        // w=0: a~r1 (ends ignored) and d~r4; b reaches r2 only at w=2 (r1 is taken by a).
        assert_eq!([0, 1, 2, 3].map(|w| tp(w, "multi_exon")), [1, 1, 2, 2]);
        assert_eq!([0, 1, 2, 3].map(|w| tp(w, "mono_exon")), [1, 1, 1, 1]);
        assert_eq!([0, 1, 2, 3].map(|w| tp(w, "")), [2, 2, 3, 3]);
        assert_eq!(results[0]["fp"], 4);
        assert_eq!(results[0]["fn"], 3);

        // Summary: same results, input and evaluated GTFs side by side, no index summaries.
        let summary: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(
            dir.join("out").join("results_summary.json"),
        )?)?;
        assert_eq!(summary["preset"], "junction");
        assert_eq!(summary["evaluated_on"], "collapsed");
        assert_eq!(
            summary["params"],
            serde_json::json!({"wob": 3, "tss_wob": null, "tes_wob": null, "mono_ovlp": 0.9})
        );
        assert_eq!(
            summary["reference"],
            serde_json::json!({
                "sample_id": "REF",
                "input_path": dir.join("ref.gtf"),
                "evaluated_path": dir.join("out/collapsed/REF/REF.merged.gtf.gz"),
                "input_total": 7, "total": 5, "multi_exon_total": 3, "mono_exon_total": 2
            })
        );
        let sq = &summary["queries"][0];
        assert_eq!(summary["query_count"], 1);
        assert_eq!(
            (&sq["sample_id"], &sq["input_total"], &sq["total"]),
            (&"Q1".into(), &7.into(), &6.into())
        );
        assert_eq!(
            (&sq["multi_exon_total"], &sq["mono_exon_total"]),
            (&4.into(), &2.into())
        );
        assert_eq!(sq["input_path"], serde_json::json!(dir.join("q1.gtf")));
        assert_eq!(sq["results"], q1["results"]);
        assert!(sq.get("gtf_stats").is_none() && sq.get("collapsed").is_none());

        // TSV: header plus one row per w, all-record metrics.
        let tsv = std::fs::read_to_string(dir.join("out").join("results_summary.tsv"))?;
        let rows: Vec<&str> = tsv.lines().collect();
        assert_eq!(rows.len(), 1 + results.len());
        assert!(rows[0].starts_with("sample_id\tinput_path\tw\t"));
        let w0: Vec<&str> = rows[1].split('\t').collect();
        assert_eq!(
            (w0[0], w0[2], w0[3], w0[5], w0[6]),
            ("Q1", "0", "2", "4", "3")
        );

        // Collapsed inputs live under collapsed/<LABEL>/ with their indexes.
        let out = dir.join("out");
        for label in ["REF", "Q1"] {
            let gtf = out
                .join("collapsed")
                .join(label)
                .join(format!("{label}.merged.gtf.gz"));
            assert!(gtf.is_file(), "{}", gtf.display());
            let mut info = gtf.clone();
            info.add_extension("index_info.json");
            assert!(info.is_file(), "{}", info.display());
        }
        // Status GTF records are collapsed records of the evaluated GTF.
        let collapsed_q1 = out.join("collapsed/Q1/Q1.merged.gtf.gz");
        let mut collapsed_ids = HashSet::new();
        for group in TxReader::new(&collapsed_q1)? {
            collapsed_ids.insert(group?.tx_id);
        }
        let read = |status: &str| read_status_gtf(&out.join("Q1"), status);
        let (tp_query, fp) = (read("tp-query")?, read("fp")?);
        assert_eq!((tp_query.len(), fp.len()), (2, 4));
        assert!(
            tp_query
                .keys()
                .chain(fp.keys())
                .all(|id| collapsed_ids.contains(id))
        );
        assert_eq!(collapsed_ids.len(), 6);
        let file = File::open(out.join("Q1").join("fp.gtf.gz"))?;
        let mut text = String::new();
        std::io::Read::read_to_string(&mut GzDecoder::new(BufReader::new(file)), &mut text)?;
        assert!(text.contains("preset=\"junction\";"));
        assert!(text.contains(&format!(
            "##ISOM <SAMPLE> id=\"Q1\"; input=\"{}\"; raw_input=\"{}\"; role=\"query\"; total=\"6\";",
            collapsed_q1.display(),
            dir.join("q1.gtf").display()
        )));

        // Collapse does not depend on --wob: same denominators.
        let w0 = run_cli(
            &dir,
            "out_w0",
            "ref.gtf",
            &["q1.gtf"],
            &["--preset", "junction"],
        )?;
        assert_eq!(w0["reference"]["total"], 5);
        assert_eq!(w0["queries"][0]["total"], 6);
        assert_eq!(w0["queries"][0]["results"][0], results[0]);

        // PLAN example: two reference records share a chain; the query equals one of them.
        write_gtf(&dir.join("ref_pair.gtf"), &REF[..2])?;
        write_gtf(&dir.join("q_one.gtf"), &REF[..1])?;
        let default = run_cli(&dir, "pair_default", "ref_pair.gtf", &["q_one.gtf"], &[])?;
        let r = &default["queries"][0]["results"][0];
        assert_eq!((&r["precision"], &r["recall"]), (&1.0.into(), &0.5.into()));
        let summary: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(
            dir.join("pair_default").join("results_summary.json"),
        )?)?;
        assert_eq!(summary["evaluated_on"], "input");
        assert_eq!(
            (&summary["params"]["tss_wob"], &summary["params"]["tes_wob"]),
            (&50.into(), &50.into())
        );
        let reference = &summary["reference"];
        assert_eq!(reference["evaluated_path"], reference["input_path"]);
        assert_eq!(
            (&reference["input_total"], &reference["total"]),
            (&2.into(), &2.into())
        );
        assert!(!dir.join("pair_default").join("collapsed").exists());
        let junction = run_cli(
            &dir,
            "pair_junction",
            "ref_pair.gtf",
            &["q_one.gtf"],
            &["--preset", "junction"],
        )?;
        let r = &junction["queries"][0]["results"][0];
        assert_eq!(
            (&r["precision"], &r["recall"], &r["f1"]),
            (&1.0.into(), &1.0.into(), &1.0.into())
        );
        std::fs::remove_dir_all(dir)?;
        Ok(())
    }

    #[test]
    fn tx_reader_restores_sources_and_reports_errors() -> AnyResult<()> {
        let path = std::env::temp_dir().join(format!(
            "isomatch-tx-reader-{}-{}.gtf.gz",
            std::process::id(),
            SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
        ));
        let open = |gtf: &str| -> Result<TxReader, BenchError> {
            let mut writer = GzEncoder::new(File::create(&path)?, Compression::default());
            writer.write_all(gtf.as_bytes())?;
            writer.finish()?;
            TxReader::new(&path)
        };
        let gtf = concat!(
            "# merged GTF\n\n",
            "chr1\tisomatch\ttranscript\t100\t600\t.\t+\t.\ttranscript_id \"m1\"; ISOM_EXONS \"3\"; ISOM_SRC \"S1:ref1:g1:100:600:ALLC:0:0:no_diff|S1:ref2:g1:100:600:ALLC:0:0:no_diff|S2:query:g2:110:590:NOTC:0:0:(1,-2,3),(2,4,-5)\";\n",
            "chr1\tisomatch\texon\t100\t200\t.\t+\t.\ttranscript_id \"m1\";\n",
            "chr1\tisomatch\texon\t300\t400\t.\t+\t.\ttranscript_id \"m1\";\n",
            "chr1\tisomatch\texon\t500\t600\t.\t+\t.\ttranscript_id \"m1\";\n",
            "chr2\tisomatch\ttranscript\t700\t800\t.\t-\t.\ttranscript_id \"m2\"; ISOM_EXONS \"1\"; ISOM_SRC \"S3:minus:g3:710:790:MONO:0:0:no_diff\";\n",
            "chr2\tisomatch\texon\t700\t800\t.\t-\t.\ttranscript_id \"m2\";"
        );
        let header = concat!(
            "##ISOM <VERSION> version=\"test\";\r\n",
            "##ISOM <SAMPLE> id=\"S1\"; input=\"/data/ref.gtf.gz\";\n",
            "\n",
            "##ISOM <SAMPLE> id=\"S3\"; input=\"minus.gtf.gz\";\n",
            r#"##ISOM <SAMPLE> id="S2"; input="/data/query \"quoted\" \\name.gtf.gz";"#,
            "\n",
            "##ISOM <FORMAT> ID=\"ISOM_SRC\"; Description=\"sources\";\n",
        );
        let mut reader = open(&format!("{header}{gtf}"))?;
        assert_eq!(
            reader.headerview,
            format!("{}# merged GTF\n", header.replace("\n\n", "\n"))
        );
        assert_eq!(
            reader.sample_list,
            vec![
                (1, "ref.gtf.gz".to_string()),
                (3, "minus.gtf.gz".to_string()),
                (2, r#"query "quoted" \name.gtf.gz"#.to_string()),
            ]
        );
        assert_eq!(reader.line_no, 8);
        let group = reader.next().unwrap()?;
        assert_eq!((&*group.chr, &*group.tx_id), ("chr1", "m1"));
        assert_eq!(group.strand, Strand::Plus);
        assert_eq!(group.ref_tx_cnt, 2);
        assert_eq!(group.n_exons, 3);
        assert_eq!(group.tx_structures.len(), 3);
        assert_eq!(group.tx_structures[0].sample_idx, 0);
        assert_eq!(group.tx_structures[0].junctions, [200, 300, 400, 500]);
        let query = &group.tx_structures[2];
        assert_eq!((query.sample_idx, query.tss, query.tes), (1, 110, 590));
        assert_eq!((&*query.orig_txid, &*query.orig_geneid), ("query", "g2"));
        assert_eq!(query.junctions, [202, 297, 396, 505]);
        assert_eq!(query.exons(), [(110, 202), (297, 396), (505, 590)]);
        let group = reader.next().unwrap()?;
        assert_eq!(group.strand, Strand::Minus);
        assert_eq!(group.ref_tx_cnt, 0);
        assert_eq!(group.n_exons, 1);
        assert_eq!(group.tx_structures.len(), 1);
        let minus = &group.tx_structures[0];
        assert_eq!((minus.sample_idx, minus.tss, minus.tes), (2, 790, 710));
        assert!(minus.junctions.is_empty());
        assert_eq!(minus.exons(), [(710, 790)]);
        assert!(reader.next().is_none());
        assert!(reader.next().is_none());
        drop(reader);

        let mut reader = open("# empty GTF")?;
        assert_eq!(reader.headerview, "# empty GTF");
        assert!(reader.sample_list.is_empty());
        assert!(reader.next().is_none());
        drop(reader);
        let mut reader = open("")?;
        assert!(reader.headerview.is_empty());
        assert!(reader.sample_list.is_empty());
        assert!(reader.next().is_none());
        drop(reader);
        let mut reader = open(gtf.trim_start_matches("# merged GTF\n\n"))?;
        assert!(reader.headerview.is_empty());
        assert!(reader.sample_list.is_empty());
        assert_eq!(reader.next().unwrap()?.ref_tx_cnt, 2);
        drop(reader);
        // Unescaped IDs containing ':' keep the transcript ID intact.
        let mut reader = open(&gtf.replace("S2:query:g2:", "S2:q:1:g2:"))?;
        let query = &reader.next().unwrap()?.tx_structures[2];
        assert_eq!((&*query.orig_txid, &*query.orig_geneid), ("q:1", "g2"));
        drop(reader);
        for invalid_header in [
            "##ISOM <SAMPLE> input=\"ref.gtf\";\n",
            "##ISOM <SAMPLE> id=\"S1\";\n",
            "##ISOM <SAMPLE> id=\"S0\"; input=\"ref.gtf\";\n",
            "##ISOM <SAMPLE> id=\"1\"; input=\"ref.gtf\";\n",
            "##ISOM <SAMPLE> id=\"Sbad\"; input=\"ref.gtf\";\n",
            "##ISOM <SAMPLE> id=\"S1\"; input=\"unterminated\n",
        ] {
            assert!(matches!(
                open(invalid_header),
                Err(BenchError::MalformedGTF { .. })
            ));
        }
        for invalid in [
            gtf.replace("ISOM_SRC", "MISSING_SRC"),
            gtf.replace("S1:", "S0:"),
            gtf.replace("ISOM_EXONS \"3\"", "ISOM_EXONS \"invalid\""),
            gtf.replace("(1,-2,3)", "(3,-2,3)"),
            gtf.replace("(1,-2,3)", "(1,invalid,3)"),
            gtf.replace("(1,-2,3)", "(1,300,3)"),
            gtf.replace("ISOM_EXONS \"3\"", "ISOM_EXONS \"4\""),
            gtf.lines().take(4).collect::<Vec<_>>().join("\n"),
            "malformed GTF\n".to_string(),
        ] {
            let mut reader = open(&invalid)?;
            let err = reader.next().unwrap().err().unwrap();
            match err {
                BenchError::ReadGTF { line_no, source } => {
                    assert!(line_no > 0);
                    assert!(matches!(
                        *source,
                        BenchError::MalformedGTF { .. }
                            | BenchError::ParseIntError(_)
                            | BenchError::LibGtf(_)
                            | BenchError::Tool(_)
                    ));
                }
                err => panic!("Expected ReadGTF, got {err:?}"),
            }
            assert!(reader.next().is_none());
        }
        std::fs::write(&path, b"invalid gzip")?;
        assert!(matches!(TxReader::new(&path), Err(BenchError::Io(_))));
        std::fs::remove_file(&path)?;
        assert!(matches!(
            TxReader::new(&path),
            Err(BenchError::OpenGTF { path: error_path, source })
                if error_path == path && source.kind() == std::io::ErrorKind::NotFound
        ));
        Ok(())
    }
}
