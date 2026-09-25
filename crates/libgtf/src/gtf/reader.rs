use std::{
    cmp::Ordering,
    cmp::Reverse,
    collections::BinaryHeap,
    fs::{self},
    io::{self, BufRead},
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

use crate::{
    error::Error,
    gtf::{
        Strand, attribute,
        bucket::{Bucket, TmpTxRec},
    },
    io::bytes_to_string,
};
use crate::{error::Result, gtf::parse_record, io::open_file_bufread};

use crate::gtf::bucket::*;
use log::warn;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::gtf::bucket::SortedBucket;

pub type ChromID = u16;

const BUCKET_COUNT: usize = 256;

/// gtf characteristics
#[derive(Debug, Clone)]
pub struct GtfProfile {
    pub chrom_names: Vec<String>,
    pub chrom_name_to_id: FxHashMap<String, ChromID>,
    pub md5: [u8; 16],
    pub file_size: u64,
}

/// GTF tx record.
/// contains everything for isomx and isoms
#[derive(Debug, Clone)]
pub struct Transcript {
    pub gidx: u64,
    pub chrom_id: ChromID,
    pub start: u32,
    pub end: u32,
    pub strand: Strand,
    pub exons: Vec<(u32, u32)>,
    pub tx_id: String,
    pub gene_id: String,
    pub is_empty: bool,
    pub attr_string: Option<Vec<u8>>,
}

impl Default for Transcript {
    fn default() -> Self {
        Self {
            gidx: 0,
            chrom_id: 0,
            start: 0,
            end: 0,
            strand: Strand::Unknown,
            exons: Vec::new(),
            tx_id: String::new(),
            gene_id: String::new(),
            is_empty: true,
            attr_string: None,
        }
    }
}

impl Transcript {
    pub fn set_gidx(&mut self, idx: u64) {
        self.gidx = idx;
    }

    pub fn set_chrom_id(&mut self, chrom_id: ChromID) {
        self.chrom_id = chrom_id;
        self.is_empty = false;
    }

    pub fn set_start(&mut self, start: u32) {
        self.start = start;
        self.is_empty = false;
    }

    pub fn get_raw_start(&self) -> u32 {
        self.start
    }

    pub fn get_0based_start(&self) -> u32 {
        self.start - 1
    }

    pub fn set_end(&mut self, end: u32) {
        self.end = end;
        self.is_empty = false;
    }

    pub fn set_strand(&mut self, strand: Strand) {
        self.strand = strand;
        self.is_empty = false;
    }

    pub fn set_tx_id(&mut self, tx_id: String) {
        self.tx_id = tx_id;
        self.is_empty = false;
    }

    pub fn set_gene_id(&mut self, gene_id: String) {
        self.gene_id = gene_id;
        self.is_empty = false;
    }

    pub fn add_exon(&mut self, exon: (u32, u32)) {
        if exon.0 < self.start || self.is_empty {
            self.start = exon.0;
        }
        if exon.1 > self.end {
            self.end = exon.1;
        }
        self.exons.push(exon);
        self.is_empty = false;
    }

    pub fn sort_exons(&mut self) {
        self.exons.sort_by_key(|e| e.0);
    }

    /// Return 0-based exon offsets relative to the left-most exon start.
    pub fn get_0based_exon_relative_offset(&self) -> Vec<(u32, u32)> {
        let base = self.exons[0].0;
        self.exons
            .iter()
            .map(|item| (item.0 - base, item.1 - base + 1))
            .collect()
    }
}

pub struct GtfReader {
    sorted_buckets: Vec<SortedBucket>,
    heap: BinaryHeap<Reverse<HeapItem>>,
    profile: GtfProfile,
    tx_counts_by_chrom_id: Vec<u64>,
    temp_dir: PathBuf,
}

impl GtfReader {
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let temp_dir = make_temp_dir(path)?;
        Self::new_with_temp_dir(path, temp_dir)
    }

    pub fn new_in<P: AsRef<Path>, Q: AsRef<Path>>(input: P, temp_parent: Q) -> Result<Self> {
        let path = input.as_ref();
        let temp_dir = make_temp_dir_in(temp_parent.as_ref())?;
        Self::new_with_temp_dir(path, temp_dir)
    }

    fn new_with_temp_dir(path: &Path, temp_dir: PathBuf) -> Result<Self> {
        let file_size = fs::metadata(path)?.len();
        let mut temp_cleanup = TempDirCleanup::new(temp_dir.clone());
        let bucket_dir = temp_dir.join("buckets");
        let sorted_dir = temp_dir.join("sorted_buckets");
        fs::create_dir(&bucket_dir)?;
        fs::create_dir(&sorted_dir)?;

        let bucket_paths = numbered_paths(&bucket_dir, "bucket", "events", BUCKET_COUNT);
        let sorted_paths = numbered_paths(&sorted_dir, "bucket", "sorted", BUCKET_COUNT);

        let mut buckets = bucket_paths
            .iter()
            .map(Bucket::init_writer)
            .collect::<io::Result<Vec<_>>>()?;

        let mut bufreader = open_file_bufread(path)?;
        let mut hasher = xxhash_rust::xxh3::Xxh3::new();
        let mut transcript_chroms: FxHashSet<String> = FxHashSet::default();
        let mut exon_chroms: FxHashSet<String> = FxHashSet::default();
        let mut has_transcript = false;
        let mut line = String::new();
        let mut line_no = 0usize;

        loop {
            line.clear();
            if bufreader.read_line(&mut line)? == 0 {
                break;
            }
            line_no += 1;
            hasher.update(line.as_bytes());

            if line.starts_with('#') {
                continue;
            }

            let line_trimmed = line.trim_end_matches(['\r', '\n']);

            let gtf_record = parse_record(line_trimmed).map_err(|err| Error::InvalidGtfLine {
                line: line_no,
                reason: err.to_string(),
            })?;
            let tx_id = attribute(gtf_record.attributes, "transcript_id").unwrap_or_default();
            let gene_id = attribute(gtf_record.attributes, "gene_id").unwrap_or_default();

            if gtf_record.feature != "transcript" && gtf_record.feature != "exon" {
                continue;
            }

            if tx_id.is_empty() || gene_id.is_empty() {
                let missing = match (tx_id.is_empty(), gene_id.is_empty()) {
                    (true, true) => "transcript_id and gene_id",
                    (true, false) => "transcript_id",
                    (false, true) => "gene_id",
                    (false, false) => unreachable!(),
                };
                return Err(Error::InvalidGtfLine {
                    line: line_no,
                    reason: format!(
                        "Missing required GTF attribute(s): {missing}. Affected line: {}",
                        line_trimmed
                    ),
                });
            }

            let hash = xxhash_rust::xxh3::xxh3_64(tx_id.as_bytes());
            let bucket_idx = (hash as usize) % BUCKET_COUNT;

            match gtf_record.feature {
                "transcript" => {
                    has_transcript = true;
                    transcript_chroms.insert(gtf_record.seqid.to_string());
                    let attr = raw_attr_bytes(line_trimmed);
                    buckets[bucket_idx].write_tx(&TmpTxRec {
                        hash,
                        chrom: gtf_record.seqid.to_owned().into_bytes(),
                        start: gtf_record.start,
                        end: gtf_record.end,
                        strand: gtf_record.strand,
                        tx_id: tx_id.to_owned().into_bytes(),
                        gene_id: gene_id.to_owned().into_bytes(),
                        raw_attr_string: attr,
                    })?;
                }
                "exon" => {
                    exon_chroms.insert(gtf_record.seqid.to_string());
                    if gtf_record.start > gtf_record.end {
                        warn!(
                            "Invalid GTF record with start > end, affected line: {}",
                            line_trimmed
                        );
                        continue;
                    }
                    buckets[bucket_idx].write_exon(&TmpExonRec {
                        hash,
                        chrom: gtf_record.seqid.to_owned().into_bytes(),
                        start: gtf_record.start,
                        end: gtf_record.end,
                        strand: gtf_record.strand,
                        tx_id: tx_id.to_owned().into_bytes(),
                        gene_id: gene_id.to_owned().into_bytes(),
                    })?;
                }
                _ => unreachable!(),
            }
        }

        // flush and close bucket writers before generate sorted buckets
        for bucket in &mut buckets {
            bucket.flush()?;
        }
        drop(buckets);

        if !has_transcript {
            return Err(Error::MissingTranscriptRecord);
        }

        //  make sure that exons and transcripts are on the same set of chrs.
        if transcript_chroms != exon_chroms {
            let mut transcript_only: Vec<String> = transcript_chroms
                .difference(&exon_chroms)
                .cloned()
                .collect();
            let mut exon_only: Vec<String> = exon_chroms
                .difference(&transcript_chroms)
                .cloned()
                .collect();
            transcript_only.sort();
            exon_only.sort();
            return Err(crate::error::Error::TranscriptExonChromMismatch {
                transcript_only,
                exon_only,
            });
        }

        // recode the chrnames to id projection
        let mut chrom_names: Vec<String> = transcript_chroms.into_iter().collect();
        chrom_names.sort();
        let chrom_name_to_id: FxHashMap<String, ChromID> = chrom_names
            .iter()
            .enumerate()
            .map(|(idx, chrom)| (chrom.clone(), (idx + 1) as ChromID))
            .collect();
        let profile = GtfProfile {
            chrom_names,
            chrom_name_to_id,
            md5: hasher.digest128().to_le_bytes(),
            file_size,
        };

        let mut tx_counts_by_chrom_id = vec![0u64; profile.chrom_names.len() + 1];
        for (bucket_path, sorted_path) in bucket_paths.iter().zip(sorted_paths.iter()) {
            let mut txs = aggregate_bucket(bucket_path, &profile)?;
            txs.sort_by(tx_sort_cmp);

            let mut sorted_bucket = SortedBucket::init_writer(sorted_path)?;
            for tx in &txs {
                if let Some(count) = tx_counts_by_chrom_id.get_mut(tx.chrom_id as usize) {
                    *count += 1;
                }
                sorted_bucket.dump_tx_structure(tx)?;
            }
            sorted_bucket.flush()?;
        }

        for path in &bucket_paths {
            let _ = fs::remove_file(path);
        }

        let mut sorted_buckets = sorted_paths
            .iter()
            .map(SortedBucket::init_reader)
            .collect::<io::Result<Vec<_>>>()?;
        let mut heap = BinaryHeap::new();
        for (bucket_idx, bucket) in sorted_buckets.iter_mut().enumerate() {
            if let Some(tx) = bucket.read_one()? {
                heap.push(Reverse(HeapItem { tx, bucket_idx }));
            }
        }

        temp_cleanup.disarm();
        Ok(Self {
            sorted_buckets,
            heap,
            profile,
            tx_counts_by_chrom_id,
            temp_dir,
        })
    }

    pub fn profile(&self) -> &GtfProfile {
        &self.profile
    }

    pub fn chrom_name(&self, chrom_id: ChromID) -> Option<&str> {
        chrom_id
            .checked_sub(1)
            .and_then(|idx| self.profile.chrom_names.get(idx as usize))
            .map(String::as_str)
    }

    pub fn transcript_count_excluding(
        &self,
        skipped_chrom_ids: &std::collections::HashSet<ChromID>,
    ) -> u64 {
        self.tx_counts_by_chrom_id
            .iter()
            .enumerate()
            .filter(|(idx, _)| !skipped_chrom_ids.contains(&(*idx as ChromID)))
            .map(|(_, count)| count)
            .sum()
    }

    /// k-way merge based on all sorted buckets
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> Result<Option<Transcript>> {
        let Some(Reverse(item)) = self.heap.pop() else {
            return Ok(None);
        };

        if let Some(next_tx) = self.sorted_buckets[item.bucket_idx].read_one()? {
            self.heap.push(Reverse(HeapItem {
                tx: next_tx,
                bucket_idx: item.bucket_idx,
            }));
        }

        Ok(Some(item.tx))
    }
}

impl Drop for GtfReader {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.temp_dir);
    }
}

struct HeapItem {
    tx: Transcript,
    bucket_idx: usize,
}

impl Eq for HeapItem {}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        tx_sort_cmp(&self.tx, &other.tx) == Ordering::Equal && self.bucket_idx == other.bucket_idx
    }
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        tx_sort_cmp(&self.tx, &other.tx).then_with(|| self.bucket_idx.cmp(&other.bucket_idx))
    }
}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn aggregate_bucket(path: &Path, profile: &GtfProfile) -> Result<Vec<Transcript>> {
    let mut bucket = Bucket::init_reader(path)?;
    let mut txs: FxHashMap<String, Rec2TxStrctureTmp> = FxHashMap::default();

    while let Some(record) = bucket.read_one()? {
        match record {
            TmpRec::Tx(tx) => observe_tmp_tx(tx, profile, &mut txs)?,
            TmpRec::Exon(exon) => observe_tmp_exon(exon, profile, &mut txs)?,
        }
    }

    let mut ready = Vec::new();
    for mut acc in txs.into_values() {
        let Some(mut tx) = acc.tx.take() else {
            continue;
        };
        if let Some(record_chrom_id) = acc.transcript_chrom_id
            && record_chrom_id != tx.chrom_id
        {
            return Err(Error::InconsistentTranscript {
                tx_id: tx.tx_id.clone(),
                reason: format!(
                    "Transcript {} record is on chrom_id {}, but its exons are on chrom_id {}",
                    tx.tx_id, record_chrom_id, tx.chrom_id
                ),
            });
        }
        tx.attr_string = acc.attr_string.take();
        tx.sort_exons();
        ready.push(tx);
    }

    Ok(ready)
}

fn observe_tmp_tx(
    record: TmpTxRec,
    profile: &GtfProfile,
    txs: &mut FxHashMap<String, Rec2TxStrctureTmp>,
) -> Result<()> {
    let tx_id = bytes_to_string(record.tx_id)?;
    let chrom_id = chrom_id_for_bytes(profile, record.chrom)?;

    let acc = txs.entry(tx_id.clone()).or_default();

    // check if the acc's chr is same as this record
    if let Some(prev_chrom_id) = acc.transcript_chrom_id {
        if prev_chrom_id != chrom_id {
            return Err(Error::InconsistentTranscript {
                tx_id: tx_id.clone(),
                reason: format!(
                    "Transcript {tx_id} has transcript records on multiple chrom_ids: {prev_chrom_id} vs {chrom_id}"
                ),
            });
        }
    } else {
        acc.transcript_chrom_id = Some(chrom_id);
    }

    acc.attr_string = Some(record.raw_attr_string);
    Ok(())
}

fn observe_tmp_exon(
    record: TmpExonRec,
    profile: &GtfProfile,
    txs: &mut FxHashMap<String, Rec2TxStrctureTmp>,
) -> Result<()> {
    let tx_id = bytes_to_string(record.tx_id)?;
    let gene_id = bytes_to_string(record.gene_id)?;
    let chrom_id = chrom_id_for_bytes(profile, record.chrom)?;
    let acc = txs.entry(tx_id.clone()).or_default();

    if let Some(record_chrom_id) = acc.transcript_chrom_id
        && record_chrom_id != chrom_id
    {
        return Err(Error::InconsistentTranscript {
            tx_id: tx_id.clone(),
            reason: format!(
                "Transcript {tx_id} record is on chrom_id {record_chrom_id}, but an exon is on chrom_id {chrom_id}"
            ),
        });
    }

    let tx = acc.tx.get_or_insert_with(|| {
        let mut tx = Transcript::default();
        tx.set_start(record.start);
        tx.set_end(record.end);
        tx.set_chrom_id(chrom_id);
        tx.set_strand(record.strand);
        tx.set_tx_id(tx_id.clone());
        tx.set_gene_id(gene_id.clone());
        tx
    });

    if tx.chrom_id != chrom_id {
        return Err(Error::InvalidGtfChr {
            reason: format!(
                "Transcript {tx_id} has exons on multiple chrom_ids: {} vs {chrom_id}",
                tx.chrom_id
            ),
        });
    }

    if tx.gene_id != gene_id {
        warn!(
            "Transcript {} has inconsistent gene_id in its exon record. chrom_id {}: {} vs {}",
            tx_id, chrom_id, tx.gene_id, gene_id
        );
    }

    if tx.strand != record.strand {
        warn!(
            "Transcript {} has inconsistent strand in its exon record at chrom_id {}. {} vs {}",
            tx_id, chrom_id, tx.strand, record.strand
        );
    }

    tx.add_exon((record.start, record.end));
    Ok(())
}

fn tx_sort_cmp(a: &Transcript, b: &Transcript) -> Ordering {
    (
        a.chrom_id,
        a.start,
        a.end,
        a.strand,
        a.tx_id.as_str(),
        a.gene_id.as_str(),
    )
        .cmp(&(
            b.chrom_id,
            b.start,
            b.end,
            b.strand,
            b.tx_id.as_str(),
            b.gene_id.as_str(),
        ))
}

fn chrom_id_for_bytes(profile: &GtfProfile, chrom: Vec<u8>) -> Result<ChromID> {
    let chrom = bytes_to_string(chrom)?;
    profile
        .chrom_name_to_id
        .get(&chrom)
        .copied()
        .ok_or_else(|| {
            // Err(
            Error::InvalidGtfChr {
                reason: format!("chromosome {chrom} was not found in GTF profile"),
            }
            // )

            // io::Error::new(
            //     ErrorKind::InvalidData,
            //     format!("chromosome {chrom} was not found in GTF profile"),
            // )
        })
}

fn raw_attr_bytes(line: &str) -> Vec<u8> {
    line.splitn(9, '\t')
        .nth(8)
        .unwrap_or("")
        .as_bytes()
        .to_vec()
}

fn numbered_paths(dir: &Path, prefix: &str, ext: &str, count: usize) -> Vec<PathBuf> {
    (0..count)
        .map(|idx| dir.join(format!("{prefix}_{idx:05}.{ext}")))
        .collect()
}

fn make_temp_dir(input_path: &Path) -> io::Result<PathBuf> {
    let parent = input_path.parent().unwrap_or_else(|| Path::new("."));
    make_temp_dir_in(parent)
}

fn make_temp_dir_in(parent: &Path) -> io::Result<PathBuf> {
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|err| io::Error::other(err.to_string()))?
        .as_nanos();
    let dir = parent.join(format!(".isomatch-index-{}-{stamp}", std::process::id()));
    fs::create_dir(&dir)?;
    Ok(dir)
}
