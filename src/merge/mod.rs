use std::{collections::HashSet, fs::File};

use crate::core::ptir::PTIR;
use crate::core::tx_base::TxBase;
use crate::core::tx_strand::ISOMSTRAND;
use crate::index::reader::ChromBlockReader;
use crate::merge::policy::merge_cluster;
use crate::{MergeArgs, index::reader::IndexReader, traits::ArgValidate};
pub mod merge_error;
pub mod mptir;
pub mod policy;
use anyhow::Context;
use anyhow::Result as AnyResult;
use anyhow::anyhow;
use log::info;
use merge_error::MergeError;
use rustc_hash::FxHashMap;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

impl ArgValidate for MergeArgs {
    fn validate(&self) {
        // place holder
    }
}

pub fn run_merge(args: MergeArgs) -> AnyResult<()> {
    // open all files (isomx) into a vec
    args.validate();
    let n_inputs = args.inputs.len();
    info!("Loading {n_inputs} gtf(s)...");
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

    info!("Loading chromosome names...");
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

    // init min-heap

    // get chromsome names and get chromblockreader from all indexxreader, for each chromsome do:

    for chrom_name in &chrom_names {
        info!("Merging chromosome {}", chrom_name);
        let mut chrom_block_readers: Vec<ChromBlockReader> = Vec::new();
        for reader in &mut fhs {
            if let std::result::Result::Ok(chrom_block_reader) =
                reader.get_chromosome_reader(chrom_name)
            {
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

            process_super_cluster(&mut super_cluster, &args);
            super_cluster.clear();
            cluster_max_end = ptir.end;
            super_cluster.push(ptir); // first ptir for next super cluster
        }
        process_super_cluster(&mut super_cluster, &args);
        // process the last super cluster

        // build the junction-level custer

        // merge canonical tx

        // non-canoniical tx to canonical tx

        // report to unified GTF
    }
    info!("Fnished!");
    Ok(())
}

pub fn process_super_cluster(super_cluster: &mut Vec<PTIR>, args: &MergeArgs) {
    println!("super cluster size {}", super_cluster.len());

    // build junc cluster
    // cluter has same strand and junction number, which is the merge unit
    let mut clusters: std::collections::HashMap<
        (ISOMSTRAND, u16),
        Vec<usize>,
        rustc_hash::FxBuildHasher,
    > = FxHashMap::default();
    for (ptir_idx, ptir) in super_cluster.iter().enumerate() {
        // println!(
        //     "{}<{},{},n={}>",
        //     ptir.tx_boundary,
        //     ptir.source_file_id,
        //     ptir.source_geneid,
        //     ptir.n_exons,
        // );
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
        // println!("");
        for sclu_idx in sclu_idxs {
            println!(
                "junc_cluster {}|{}, ptir: {}",
                n_exons, strand, super_cluster[*sclu_idx]
            );
        }

        merge_cluster(*n_exons, *strand, sclu_idxs, super_cluster, args);
    }
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
            if let Some(tx_base) = reader.next()? {
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
            (state.next()?, state.file_id)
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
