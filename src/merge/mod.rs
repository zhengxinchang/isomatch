use crate::IndexArgs;
use crate::constants::ISOM_GTF_SCHEMA;
use crate::core::ptir::PTIR;
use crate::core::tx_base::TxBase;
use crate::core::tx_strand::ISOMSTRAND;
use crate::index::reader::ChromBlockReader;
use crate::index::run_index;
use crate::merge::grouped_ptirs::GroupedPTIR;
use crate::merge::guide::GuideDb;
use crate::merge::policy::MergePolicyUsed;
use crate::merge::policy::merge_cluster;
use crate::utils::check_index_ready;
use crate::utils::greetings2;
use crate::utils::print_json_block;
use crate::utils::require_file;
use crate::{MergeArgs, index::reader::IndexReader, traits::ArgValidate};
use serde::Serialize;
use std::io::BufWriter;
use std::io::Write;
use std::path::PathBuf;
use std::{
    collections::{BTreeMap, HashSet},
    fs::File,
};
pub mod grouped_ptirs;
pub mod guide;
pub mod merge_error;
pub mod policy;
use anyhow::Context;
use anyhow::Result as AnyResult;
use anyhow::anyhow;
use flate2::Compression;
use flate2::write::GzEncoder;
use log::{error, info};
use merge_error::MergeError;
use rustc_hash::FxHashMap;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

impl ArgValidate for MergeArgs {
    fn validate(&self) {
        let mut error_msg = String::new();
        let mut has_error = false;

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

        for input in &self.inputs {
            require_file("Input GTF", input, &mut error_msg, &mut has_error);
        }

        if let Some(guide_tss) = &self.guide_tss {
            require_file("TSS BED", guide_tss, &mut error_msg, &mut has_error);
        }

        if let Some(guide_tes) = &self.guide_tes {
            require_file("TES BED", guide_tes, &mut error_msg, &mut has_error);
        }

        if let Some(chrmap) = &self.chrmap {
            require_file("Chromosome map", chrmap, &mut error_msg, &mut has_error);
        }

        if has_error {
            error!("Error validating arguments: {}", error_msg);
            std::process::exit(1);
        }
    }
}

pub fn run_merge(args: MergeArgs) -> AnyResult<()> {
    greetings2(&args);
    // open all files (isomx) into a vec
    args.validate();
    let n_inputs = args.inputs.len();
    let input_file_names = input_file_name_bytes(&args.inputs);
    info!("Loading {n_inputs} gtf(s)");
    // auto indexing
    info!(
        "Checking input GTF indexes; missing/corrupted/outdated indexes will be created automatically."
    );

    for gtf in &args.inputs {
        if !check_index_ready(gtf) {
            info!("Indexing {}", gtf.display());
            let mut index_args = IndexArgs {
                input: gtf.clone(),
                ref_fa: args.ref_fa.clone(),
                seqfa: None,
                out: None,
                skip_missing_ref_chr: args.skip_missing_ref_chr,
                quiet: true,
            };
            run_index(&mut index_args)
                .with_context(|| format!("Can not index GTF: {}", gtf.display()))?;
        }
    }

    let mut stats = MergeStats {
        source_files: n_inputs as u32,
        ..MergeStats::default()
    };

    let mut gtf_out_path = args.out.clone();
    gtf_out_path.add_extension("merged.gtf.gz");

    let mut merge_info_path = args.out.clone();
    merge_info_path.add_extension("merged_info.json");

    let mut fhs: Vec<IndexReader> = Vec::with_capacity(n_inputs);
    for (file_id, input_path) in args.inputs.iter().enumerate() {
        let mut index_path = input_path.clone();
        index_path.add_extension("isomx");

        let f = File::open(&index_path)
            .with_context(|| format!("Can not open index: {}", index_path.display()))?;

        let reader = match IndexReader::open(f, file_id) {
            Ok(reader) => reader,
            Err(e) => {
                return Err(anyhow!(
                    "Can not load index {}: {}",
                    index_path.display(),
                    e
                ));
            }
        };

        fhs.push(reader);
    }

    info!("Loading chromosome names");
    // collect all chromsome from all files and build a unique list
    let mut chrom_names = Vec::new();
    let mut seen_chroms = HashSet::new();
    for reader in &fhs {
        for chrom_name in &reader.chrom_names {
            if seen_chroms.insert(chrom_name.clone()) {
                chrom_names.push(chrom_name.clone());
            }
        }
    }

    // load the guide files if provided

    let guide_tss_index = if let Some(tss_p) = args.guide_tss.clone() {
        Some(GuideDb::from_bed_path(
            tss_p,
            guide::GuideBEDType::Tss,
            &args.chrmap,
        )?)
    } else {
        None
    };

    let guide_tes_index = if let Some(tes_p) = args.guide_tes.clone() {
        Some(GuideDb::from_bed_path(
            tes_p,
            guide::GuideBEDType::Tes,
            &args.chrmap,
        )?)
    } else {
        None
    };

    let f = File::create(&gtf_out_path)?;
    let mut bufwriter: Box<dyn Write> =
        Box::new(BufWriter::new(GzEncoder::new(f, Compression::default())));

    let track_out_path = PathBuf::from(format!("{}.track.tsv.gz", args.out.display()));
    let track_f = File::create(&track_out_path)?;
    let mut track_bufwriter: Box<dyn Write> = Box::new(BufWriter::new(GzEncoder::new(
        track_f,
        Compression::default(),
    )));
    writeln!(
        track_bufwriter,
        "merged_tx_id\tmerged_gene_id\tmerged_start\tmerged_end\tmerged_strand\tmerged_exon_num\tjunction_policy\ttss_policy\ttes_policy\tsrc_tx_count_in_merged_group\tsrc_tx_id\tsrc_gene_id\ttotal_donor_diff\ttotal_acceptor_diff\texon_diff\tsrc_file_name"
    )?;

    add_output_header(&mut bufwriter, &args)?;

    // get chromsome names and get chromblockreader from all indexxreader, for each chromsome do:
    let mut global_scluster_id = 0u32;
    let mut global_tx_id = 0u32;
    for chrom_name in &chrom_names {
        info!("Merging chromosome {}", chrom_name);
        let mut chrom_block_readers: Vec<ChromBlockReader> = Vec::new();
        for reader in &mut fhs {
            if let Ok(chrom_block_reader) = reader.get_chromosome_reader(chrom_name) {
                chrom_block_readers.push(chrom_block_reader);
            }
        }

        // k-way merge, build super cluster

        let mut kway_merger = KwayMerger::new(chrom_block_readers)?;

        let mut super_cluster: Vec<PTIR> = Vec::new();

        let first_ptir = match kway_merger.try_next()? {
            Some(ptir) => ptir,
            None => continue,
        };

        let mut cluster_max_end = first_ptir.end;
        super_cluster.push(first_ptir);

        for ptir in kway_merger {
            if ptir.start <= cluster_max_end {
                cluster_max_end = cluster_max_end.max(ptir.end);
                super_cluster.push(ptir);
                continue;
            }
            // global_scluster_id += 1;
            process_super_cluster(
                chrom_name,
                &mut super_cluster,
                &mut global_scluster_id,
                &mut global_tx_id,
                &mut stats,
                &args,
                &input_file_names,
                bufwriter.as_mut(),
                track_bufwriter.as_mut(),
                &guide_tss_index,
                &guide_tes_index,
            )?;
            super_cluster.clear();
            cluster_max_end = ptir.end;
            super_cluster.push(ptir); // first ptir for next super cluster
        }
        // process the last super cluster
        // global_scluster_id += 1;
        process_super_cluster(
            chrom_name,
            &mut super_cluster,
            &mut global_scluster_id,
            &mut global_tx_id,
            &mut stats,
            &args,
            &input_file_names,
            bufwriter.as_mut(),
            track_bufwriter.as_mut(),
            &guide_tss_index,
            &guide_tes_index,
        )?;

        // report to unified GTF
    }

    bufwriter.flush()?;
    drop(bufwriter);
    track_bufwriter.flush()?;
    drop(track_bufwriter);
    stats.finalize();
    print_json_block("Merge summary", &stats);
    info!("Output file saved at: {}", gtf_out_path.to_string_lossy());
    info!("Track file saved at: {}", track_out_path.to_string_lossy());

    let mut merge_info_writer = File::create(&merge_info_path)?;

    let msg = serde_json::to_string_pretty(&stats)?;

    merge_info_writer.write(msg.as_bytes())?;
    merge_info_writer.flush()?;

    info!("Finished!");
    Ok(())
}

pub fn process_super_cluster(
    chrom_name: &str,
    super_cluster: &mut Vec<PTIR>,
    global_scluster_id: &mut u32,
    global_tx_id: &mut u32,
    stats: &mut MergeStats,
    args: &MergeArgs,
    input_file_names: &[Vec<u8>],
    bufwriter: &mut dyn Write,
    track_bufwriter: &mut dyn Write,
    guide_tss: &Option<GuideDb>,
    guide_tes: &Option<GuideDb>,
) -> Result<(), MergeError> {
    // println!("super cluster size {}", super_cluster.len());
    *global_scluster_id += 1;
    stats.observe_source_txs(super_cluster.len());
    // build junc cluster
    // cluter has same strand and junction number, which is the merge unit
    let mut clusters: std::collections::HashMap<
        (ISOMSTRAND, u16),
        Vec<usize>,
        rustc_hash::FxBuildHasher,
    > = FxHashMap::default();
    for (ptir_idx, ptir) in super_cluster.iter().enumerate() {
        let key: (ISOMSTRAND, u16) = (ptir.strand, ptir.n_exons);
        let cluster = clusters.entry(key).or_insert(Vec::new());
        cluster.push(ptir_idx);
    }

    // make sure clusters are sorted by strand and then n_exon(desending).
    // then when a cluster has small exon, it can get next
    let mut cluster_items: Vec<_> = clusters.iter().collect();
    cluster_items.sort_by_key(|((strand, n_exons), _)| (*strand, std::cmp::Reverse(*n_exons)));

    // process each junc cluster
    for ((strand, n_exons), sclu_idxs) in cluster_items {
        let mut grpptirs = merge_cluster(
            chrom_name,
            *n_exons,
            *strand,
            sclu_idxs,
            super_cluster,
            args,
            guide_tss,
            guide_tes,
        )?;

        for grpptir in grpptirs.iter_mut() {
            *global_tx_id += 1;
            grpptir.update_ids(*global_scluster_id, *global_tx_id);
            stats.observe_merged_tx(grpptir);
            grpptir.write_gtf_block(
                chrom_name,
                super_cluster,
                bufwriter,
                track_bufwriter,
                input_file_names,
            )?;
        }
    }

    Ok(())
}

pub struct KwayMerger {
    readers: Vec<ChromBlockReader>,
    heap: BinaryHeap<Reverse<(TxBase, usize, usize)>>,
}

impl KwayMerger {
    fn new(mut readers: Vec<ChromBlockReader>) -> Result<Self, MergeError> {
        let mut heap: BinaryHeap<Reverse<(TxBase, usize, usize)>> = BinaryHeap::new();

        for (idx, reader) in readers.iter_mut().enumerate() {
            let file_id = reader.file_id;
            if let Some(tx_base) = reader.next_record()? {
                heap.push(Reverse((tx_base, file_id, idx)));
            }
        }

        Ok(Self { readers, heap })
    }

    fn try_next(&mut self) -> Result<Option<PTIR>, MergeError> {
        let Some(Reverse((tx_base, _file_id, vec_idx))) = self.heap.pop() else {
            return Ok(None);
        };

        let (next_tx_base, file_id) = {
            let state = &mut self.readers[vec_idx];
            (state.next_record()?, state.file_id)
        };

        if let Some(next_tx_base) = next_tx_base {
            self.heap.push(Reverse((next_tx_base, file_id, vec_idx)));
        }

        let ptir = PTIR::from_tx_base(
            tx_base,
            file_id,
            &self.readers[vec_idx].junction_pool,
            &self.readers[vec_idx].splice_site_pool,
            &self.readers[vec_idx].string_pool,
        );

        Ok(Some(ptir))
    }
}

fn input_file_name_bytes(inputs: &[PathBuf]) -> Vec<Vec<u8>> {
    inputs
        .iter()
        .map(|input| {
            let file_name = input.file_name().unwrap_or_else(|| input.as_os_str());
            os_str_bytes(file_name)
        })
        .collect()
}

#[cfg(unix)]
fn os_str_bytes(value: &std::ffi::OsStr) -> Vec<u8> {
    use std::os::unix::ffi::OsStrExt;

    value.as_bytes().to_vec()
}

#[cfg(not(unix))]
fn os_str_bytes(value: &std::ffi::OsStr) -> Vec<u8> {
    value.to_string_lossy().as_bytes().to_vec()
}

impl Iterator for KwayMerger {
    type Item = PTIR;

    fn next(&mut self) -> Option<Self::Item> {
        match self.try_next() {
            Ok(item) => item,
            Err(e) => {
                panic!("Can not read next ptir because: {}", e);
            }
        }
    }
}
pub fn add_output_header(bufwriter: &mut dyn Write, args: &MergeArgs) -> AnyResult<()> {
    // ##ISOM <VERSION> version = 1.0; link = "github.."
    // ##ISOM <SAMPLE> ID="S1"; Name="xxx.gtf.gz"
    // ##ISOM <SAMPLE> ID="S2"; Name="xxx.gtf.gz"
    // ##ISOM <FORMAT> ISOM_COUNT = ""
    // ##ISOM <FORMAT> ISOM_SRC = ""
    // ##ISOM <COMMAND> isomatch ...

    let escape = |value: &str| value.replace('\\', "\\\\").replace('"', "\\\"");

    writeln!(
        bufwriter,
        "##ISOM <VERSION> version=\"{}\"; program=\"{}\"; schema=\"{}\"",
        env!("CARGO_PKG_VERSION"),
        env!("CARGO_PKG_NAME"),
        ISOM_GTF_SCHEMA
    )?;

    for (idx, input) in args.inputs.iter().enumerate() {
        writeln!(
            bufwriter,
            "##ISOM <SAMPLE> id=\"S{}\"; input=\"{}\";",
            idx + 1,
            escape(&input.to_string_lossy())
        )?;
    }

    writeln!(
        bufwriter,
        "##ISOM <FORMAT> ID=\"ISOM_EXONS\"; Description=\"number of exons for this transcript\";"
    )?;

    writeln!(
        bufwriter,
        "##ISOM <FORMAT> ID=\"ISOM_COUNT\"; Description=\"number of source transcripts merged into this output transcript\";"
    )?;
    writeln!(
        bufwriter,
        "##ISOM <FORMAT> ID=\"ISOM_SRC\"; Description=\"vertical line separated source transcript records in the form S#:tx_id:start:end:tx_type:donor_diff:acceptor_diff:(exon_number,left_offset,right_offset),(exon_number,left_offset,right_offset)... Only exons has difference will be shown.\";"
    )?;
    writeln!(
        bufwriter,
        "##ISOM <FORMAT> ID=\"ISOM_REPR_POLICY\"; Description=\"representative selection policies recorded as SJ_POLICY:TSS_POLICY:TES_POLICY, with NA for non-applicable fields\";"
    )?;

    let command = std::env::args_os()
        .map(|arg| escape(&arg.to_string_lossy()))
        .collect::<Vec<_>>()
        .join(" ");

    writeln!(bufwriter, "##ISOM <COMMAND> cmd={}", command)?;

    Ok(())
}

#[derive(Debug, Default, Serialize)]
pub struct MergeStats {
    pub source_files: u32,
    pub total_tx_cnt: u32,
    pub merged_tx_cnt: u32,
    pub merged_multi_exons_tx_cnt: u32,
    pub merged_mono_exon_tx_cnt: u32,
    pub tss_guide_cnt: u32,
    pub tss_guide_pct: f64,
    pub tes_guide_cnt: u32,
    pub tes_guide_pct: f64,
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub merged_tx_by_source_cnt: BTreeMap<u32, u32>,
}

impl MergeStats {
    fn observe_source_txs(&mut self, tx_count: usize) {
        self.total_tx_cnt += tx_count as u32;
    }

    fn observe_merged_tx(&mut self, grpptir: &GroupedPTIR) {
        self.merged_tx_cnt += 1;

        if grpptir.n_exon() <= 1 {
            self.merged_mono_exon_tx_cnt += 1;
        } else {
            self.merged_multi_exons_tx_cnt += 1;
        }

        if matches!(grpptir.used_tss_policy(), MergePolicyUsed::Guide(_)) {
            self.tss_guide_cnt += 1;
        }

        if matches!(grpptir.used_tes_policy(), MergePolicyUsed::Guide(_)) {
            self.tes_guide_cnt += 1;
        }

        *self
            .merged_tx_by_source_cnt
            .entry(grpptir.total_count())
            .or_insert(0) += 1;
    }

    fn finalize(&mut self) {
        if self.merged_tx_cnt == 0 {
            self.tss_guide_pct = 0.0;
            self.tes_guide_pct = 0.0;
            return;
        }

        self.tss_guide_pct =
            (self.tss_guide_cnt as f64 * 100.0 / self.merged_tx_cnt as f64 * 10000.0).round()
                / 10000.0;
        self.tes_guide_pct =
            (self.tes_guide_cnt as f64 * 100.0 / self.merged_tx_cnt as f64 * 10000.0).round()
                / 10000.0;
    }
}
