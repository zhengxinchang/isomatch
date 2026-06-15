use std::{
    cmp::Ordering,
    cmp::Reverse,
    collections::BinaryHeap,
    fs::{self, File},
    io::{self, BufRead, BufReader, BufWriter, Cursor, Error, ErrorKind, Read, Write},
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

use log::{error, warn};
use rustc_hash::{FxHashMap, FxHashSet};
use thiserror::Error;

use crate::{core::tx_strand::ISOMSTRAND, utils::open_file_bufread};

pub type ChromID = u16;

const BUCKET_COUNT: usize = 256;
const REC_TX: u8 = 1;
const REC_EXON: u8 = 2;

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
pub struct TxStructure {
    pub gidx: u64,
    pub chrom_id: ChromID,
    pub start: u32,
    pub end: u32,
    pub strand: ISOMSTRAND,
    pub exons: Vec<(u32, u32)>,
    pub tx_id: String,
    pub gene_id: String,
    pub is_empty: bool,
    pub attr_string: Option<Vec<u8>>,
}

impl TxStructure {
    pub fn default() -> Self {
        Self {
            gidx: 0,
            chrom_id: 0,
            start: 0,
            end: 0,
            strand: ISOMSTRAND::Unknown,
            exons: Vec::new(),
            tx_id: String::new(),
            gene_id: String::new(),
            is_empty: true,
            attr_string: None,
        }
    }

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

    pub fn set_strand(&mut self, strand: ISOMSTRAND) {
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

pub struct Bucket {
    writer: Option<BufWriter<File>>,
    reader: Option<BufReader<File>>,
}

impl Bucket {
    pub fn init_writer<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self {
            writer: Some(BufWriter::new(File::create(path.as_ref())?)),
            reader: None,
        })
    }

    pub fn init_reader<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self {
            writer: None,
            reader: Some(BufReader::new(File::open(path.as_ref())?)),
        })
    }

    pub fn flush(&mut self) -> io::Result<()> {
        if let Some(writer) = self.writer.as_mut() {
            writer.flush()?;
        }
        Ok(())
    }

    pub fn write_tx(&mut self, tx: &TmpTxRec) -> io::Result<()> {
        let mut payload = Vec::new();
        write_u64(&mut payload, tx.hash)?;
        write_bytes(&mut payload, &tx.chrom)?;
        write_u32(&mut payload, tx.start)?;
        write_u32(&mut payload, tx.end)?;
        write_u8(&mut payload, tx.strand as u8)?;
        write_bytes(&mut payload, &tx.tx_id)?;
        write_bytes(&mut payload, &tx.gene_id)?;
        write_bytes(&mut payload, &tx.raw_attr_string)?;
        self.write_record(REC_TX, &payload)
    }

    pub fn write_exon(&mut self, exon: &TmpExonRec) -> io::Result<()> {
        let mut payload = Vec::new();
        write_u64(&mut payload, exon.hash)?;
        write_bytes(&mut payload, &exon.chrom)?;
        write_u32(&mut payload, exon.start)?;
        write_u32(&mut payload, exon.end)?;
        write_u8(&mut payload, exon.strand as u8)?;
        write_bytes(&mut payload, &exon.tx_id)?;
        write_bytes(&mut payload, &exon.gene_id)?;
        self.write_record(REC_EXON, &payload)
    }

    pub fn read_one(&mut self) -> io::Result<Option<TmpRec>> {
        let Some(reader) = self.reader.as_mut() else {
            return Err(io::Error::new(
                ErrorKind::InvalidInput,
                "bucket is not open for reading",
            ));
        };

        let Some(kind) = read_u8_opt(reader)? else {
            return Ok(None);
        };
        let payload_len = read_u32(reader)? as usize;
        let mut payload = vec![0u8; payload_len];
        reader.read_exact(&mut payload)?;
        let mut cursor = Cursor::new(payload);

        match kind {
            REC_TX => Ok(Some(TmpRec::Tx(TmpTxRec {
                hash: read_u64(&mut cursor)?,
                chrom: read_bytes(&mut cursor)?,
                start: read_u32(&mut cursor)?,
                end: read_u32(&mut cursor)?,
                strand: read_strand(&mut cursor)?,
                tx_id: read_bytes(&mut cursor)?,
                gene_id: read_bytes(&mut cursor)?,
                raw_attr_string: read_bytes(&mut cursor)?,
            }))),
            REC_EXON => Ok(Some(TmpRec::Exon(TmpExonRec {
                hash: read_u64(&mut cursor)?,
                chrom: read_bytes(&mut cursor)?,
                start: read_u32(&mut cursor)?,
                end: read_u32(&mut cursor)?,
                strand: read_strand(&mut cursor)?,
                tx_id: read_bytes(&mut cursor)?,
                gene_id: read_bytes(&mut cursor)?,
            }))),
            _ => Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("unknown bucket record kind {kind}"),
            )),
        }
    }

    fn write_record(&mut self, kind: u8, payload: &[u8]) -> io::Result<()> {
        let Some(writer) = self.writer.as_mut() else {
            return Err(io::Error::new(
                ErrorKind::InvalidInput,
                "bucket is not open for writing",
            ));
        };
        write_u8(writer, kind)?;
        write_bytes(writer, payload)
    }
}

pub enum TmpRec {
    Tx(TmpTxRec),
    Exon(TmpExonRec),
}

pub struct TmpTxRec {
    hash: u64,
    chrom: Vec<u8>,
    start: u32,
    end: u32,
    strand: ISOMSTRAND,
    tx_id: Vec<u8>,
    gene_id: Vec<u8>,
    raw_attr_string: Vec<u8>,
}

pub struct TmpExonRec {
    hash: u64,
    chrom: Vec<u8>,
    start: u32,
    end: u32,
    strand: ISOMSTRAND,
    tx_id: Vec<u8>,
    gene_id: Vec<u8>,
}

pub struct SortedBucket {
    pub writer: Option<BufWriter<File>>,
    pub reader: Option<BufReader<File>>,
}

impl SortedBucket {
    pub fn init_writer<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self {
            writer: Some(BufWriter::new(File::create(path.as_ref())?)),
            reader: None,
        })
    }

    pub fn init_reader<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self {
            writer: None,
            reader: Some(BufReader::new(File::open(path.as_ref())?)),
        })
    }

    pub fn dump_tx_structure(&mut self, tx: &TxStructure) -> io::Result<()> {
        let Some(writer) = self.writer.as_mut() else {
            return Err(io::Error::new(
                ErrorKind::InvalidInput,
                "sorted bucket is not open for writing",
            ));
        };

        let mut payload = Vec::new();
        write_u16(&mut payload, tx.chrom_id)?;
        write_u32(&mut payload, tx.start)?;
        write_u32(&mut payload, tx.end)?;
        write_u8(&mut payload, tx.strand as u8)?;
        write_bytes(&mut payload, tx.tx_id.as_bytes())?;
        write_bytes(&mut payload, tx.gene_id.as_bytes())?;
        write_u32(&mut payload, tx.exons.len() as u32)?;
        for (start, end) in &tx.exons {
            write_u32(&mut payload, *start)?;
            write_u32(&mut payload, *end)?;
        }
        match &tx.attr_string {
            Some(attr) => {
                write_u8(&mut payload, 1)?;
                write_bytes(&mut payload, attr)?;
            }
            None => write_u8(&mut payload, 0)?,
        }
        write_bytes(writer, &payload)
    }

    pub fn flush(&mut self) -> io::Result<()> {
        if let Some(writer) = self.writer.as_mut() {
            writer.flush()?;
        }
        Ok(())
    }

    pub fn read_one(&mut self) -> io::Result<Option<TxStructure>> {
        let Some(reader) = self.reader.as_mut() else {
            return Err(io::Error::new(
                ErrorKind::InvalidInput,
                "sorted bucket is not open for reading",
            ));
        };

        let Some(payload_len) = read_u32_opt(reader)? else {
            return Ok(None);
        };
        let mut payload = vec![0u8; payload_len as usize];
        reader.read_exact(&mut payload)?;
        let mut cursor = Cursor::new(payload);

        let chrom_id = read_u16(&mut cursor)?;
        let start = read_u32(&mut cursor)?;
        let end = read_u32(&mut cursor)?;
        let strand = read_strand(&mut cursor)?;
        let tx_id = bytes_to_string_io(read_bytes(&mut cursor)?)?;
        let gene_id = bytes_to_string_io(read_bytes(&mut cursor)?)?;
        let exon_count = read_u32(&mut cursor)? as usize;
        let mut exons = Vec::with_capacity(exon_count);
        for _ in 0..exon_count {
            exons.push((read_u32(&mut cursor)?, read_u32(&mut cursor)?));
        }
        let attr_string = match read_u8(&mut cursor)? {
            0 => None,
            1 => Some(read_bytes(&mut cursor)?),
            value => {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    format!("invalid attr presence flag {value}"),
                ));
            }
        };

        Ok(Some(TxStructure {
            gidx: 0,
            chrom_id,
            start,
            end,
            strand,
            exons,
            tx_id,
            gene_id,
            is_empty: false,
            attr_string,
        }))
    }
}

/// intermediate structure for aggregating exon and transcript records in one TxStrcture
#[derive(Default)]
struct Rec2TxStrctureTmp {
    tx: Option<TxStructure>,
    transcript_chrom_id: Option<ChromID>,
    attr_string: Option<Vec<u8>>,
}

pub struct MyGTFReader {
    sorted_buckets: Vec<SortedBucket>,
    heap: BinaryHeap<Reverse<HeapItem>>,
    profile: GtfProfile,
    tx_counts_by_chrom_id: Vec<u64>,
    temp_dir: PathBuf,
}

impl MyGTFReader {
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, GTFError> {
        let path = path.as_ref();
        let file_size = fs::metadata(path)?.len();
        let temp_dir = make_temp_dir(path)?;
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
            let (chrom, feat, start, end, strand, tx_id, gene_id) = process_gtf_line(line_trimmed)
                .map_err(|err| {
                    GTFError::Io(Error::new(
                        err.kind(),
                        format!("Invalid GTF at line {line_no}: {err}"),
                    ))
                })?;

            if feat.as_str() != "transcript" && feat.as_str() != "exon" {
                continue;
            }

            if tx_id.is_empty() || gene_id.is_empty() {
                let missing = match (tx_id.is_empty(), gene_id.is_empty()) {
                    (true, true) => "transcript_id and gene_id",
                    (true, false) => "transcript_id",
                    (false, true) => "gene_id",
                    (false, false) => unreachable!(),
                };
                return Err(GTFError::Io(Error::new(
                    ErrorKind::InvalidData,
                    format!(
                        "Missing required GTF attribute(s): {missing}. Affected line: {}",
                        line_trimmed
                    ),
                )));
            }

            let hash = xxhash_rust::xxh3::xxh3_64(tx_id.as_bytes());
            let bucket_idx = (hash as usize) % BUCKET_COUNT;

            match feat.as_str() {
                "transcript" => {
                    has_transcript = true;
                    transcript_chroms.insert(chrom.clone());
                    let attr = raw_attr_bytes(line_trimmed);
                    buckets[bucket_idx].write_tx(&TmpTxRec {
                        hash,
                        chrom: chrom.into_bytes(),
                        start,
                        end,
                        strand,
                        tx_id: tx_id.into_bytes(),
                        gene_id: gene_id.into_bytes(),
                        raw_attr_string: attr,
                    })?;
                }
                "exon" => {
                    exon_chroms.insert(chrom.clone());
                    if start > end {
                        warn!(
                            "Invalid GTF record with start > end, affected line: {}",
                            line_trimmed
                        );
                        continue;
                    }
                    buckets[bucket_idx].write_exon(&TmpExonRec {
                        hash,
                        chrom: chrom.into_bytes(),
                        start,
                        end,
                        strand,
                        tx_id: tx_id.into_bytes(),
                        gene_id: gene_id.into_bytes(),
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
            return Err(GTFError::MissingTranscriptRecord);
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
            return Err(GTFError::TranscriptExonChromMismatch {
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
    pub fn next(&mut self) -> Result<Option<TxStructure>, GTFError> {
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

impl Drop for MyGTFReader {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.temp_dir);
    }
}

struct HeapItem {
    tx: TxStructure,
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

fn aggregate_bucket(path: &Path, profile: &GtfProfile) -> Result<Vec<TxStructure>, GTFError> {
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
        if let Some(record_chrom_id) = acc.transcript_chrom_id {
            if record_chrom_id != tx.chrom_id {
                return Err(GTFError::Io(Error::new(
                    ErrorKind::InvalidData,
                    format!(
                        "Transcript {} record is on chrom_id {}, but its exons are on chrom_id {}",
                        tx.tx_id, record_chrom_id, tx.chrom_id
                    ),
                )));
            }
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
) -> Result<(), GTFError> {
    let tx_id = bytes_to_string(record.tx_id)?;
    let chrom_id = chrom_id_for_bytes(profile, record.chrom)?;
    
    let acc = txs.entry(tx_id.clone()).or_default();

    // check if the acc's chr is same as this record
    if let Some(prev_chrom_id) = acc.transcript_chrom_id {
        if prev_chrom_id != chrom_id {
            return Err(GTFError::Io(Error::new(
                ErrorKind::InvalidData,
                format!(
                    "Transcript {tx_id} has transcript records on multiple chrom_ids: {prev_chrom_id} vs {chrom_id}"
                ),
            )));
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
) -> Result<(), GTFError> {
    let tx_id = bytes_to_string(record.tx_id)?;
    let gene_id = bytes_to_string(record.gene_id)?;
    let chrom_id = chrom_id_for_bytes(profile, record.chrom)?;
    let acc = txs.entry(tx_id.clone()).or_default();

    if let Some(record_chrom_id) = acc.transcript_chrom_id {
        if record_chrom_id != chrom_id {
            return Err(GTFError::Io(Error::new(
                ErrorKind::InvalidData,
                format!(
                    "Transcript {tx_id} record is on chrom_id {record_chrom_id}, but an exon is on chrom_id {chrom_id}"
                ),
            )));
        }
    }

    let tx = acc.tx.get_or_insert_with(|| {
        let mut tx = TxStructure::default();
        tx.set_start(record.start);
        tx.set_end(record.end);
        tx.set_chrom_id(chrom_id);
        tx.set_strand(record.strand);
        tx.set_tx_id(tx_id.clone());
        tx.set_gene_id(gene_id.clone());
        tx
    });

    if tx.chrom_id != chrom_id {
        return Err(GTFError::Io(Error::new(
            ErrorKind::InvalidData,
            format!(
                "Transcript {tx_id} has exons on multiple chrom_ids: {} vs {chrom_id}",
                tx.chrom_id
            ),
        )));
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

fn tx_sort_cmp(a: &TxStructure, b: &TxStructure) -> Ordering {
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

fn chrom_id_for_bytes(profile: &GtfProfile, chrom: Vec<u8>) -> Result<ChromID, GTFError> {
    let chrom = bytes_to_string(chrom)?;
    profile
        .chrom_name_to_id
        .get(&chrom)
        .copied()
        .ok_or_else(|| {
            GTFError::Io(Error::new(
                ErrorKind::InvalidData,
                format!("chromosome {chrom} was not found in GTF profile"),
            ))
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
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|err| io::Error::new(ErrorKind::Other, err.to_string()))?
        .as_nanos();
    let dir = parent.join(format!(".isomatch-index-{}-{stamp}", std::process::id()));
    fs::create_dir(&dir)?;
    Ok(dir)
}

pub fn process_gtf_line(
    s: &str,
) -> Result<(String, String, u32, u32, ISOMSTRAND, String, String), Error> {
    let parts: Vec<&str> = s.split('\t').collect();

    if parts.len() < 9 {
        return Err(Error::new(
            ErrorKind::InvalidData,
            format!(
                "Invalid GTF line: fewer than 9 columns. Affected line: {}",
                s.trim_end()
            ),
        ));
    }

    let chrom = parts[0].to_string();
    let feature_type = parts[2].to_string();
    let start = parts[3].parse::<u32>().map_err(|_| {
        Error::new(
            ErrorKind::InvalidData,
            format!("Invalid start coordinate. Affected line: {}", s.trim_end()),
        )
    })?;
    let end = parts[4].parse::<u32>().map_err(|_| {
        Error::new(
            ErrorKind::InvalidData,
            format!("Invalid end coordinate. Affected line: {}", s.trim_end()),
        )
    })?;
    let strand = match parts[6] {
        "-" => ISOMSTRAND::Minus,
        "+" => ISOMSTRAND::Plus,
        "." => ISOMSTRAND::Unknown,
        _ => {
            error!("Unknown Strand for transcript: {s} ");
            std::process::exit(1);
        }
    };
    let (tx_id, gene_id) = parse_gtf_attributes(parts[8]);

    Ok((chrom, feature_type, start, end, strand, tx_id, gene_id))
}

/// Take the attributes column of a GTF line and extract one value, supporting
/// both quoted and unquoted formats.
pub(crate) fn parse_gtf_attr_value(attrs: &str, key: &str) -> Option<String> {
    let mut saw_empty_match = false;

    for attr in attrs.split(';') {
        let attr = attr.trim();
        if attr.is_empty() || !attr.starts_with(key) {
            continue;
        }

        let value = extract_attr_value(attr);
        if !value.is_empty() {
            return Some(value);
        }
        saw_empty_match = true;
    }

    saw_empty_match.then(String::new)
}

fn parse_gtf_attributes(attrs: &str) -> (String, String) {
    (
        parse_gtf_attr_value(attrs, "transcript_id").unwrap_or_default(),
        parse_gtf_attr_value(attrs, "gene_id").unwrap_or_default(),
    )
}

fn extract_attr_value(attr: &str) -> String {
    if let Some(q_start) = attr.find('"') {
        if let Some(q_len) = attr[q_start + 1..].find('"') {
            return attr[q_start + 1..q_start + 1 + q_len].to_string();
        }
    }
    attr.split_ascii_whitespace()
        .nth(1)
        .unwrap_or("")
        .to_string()
}

fn write_u8<W: Write>(writer: &mut W, value: u8) -> io::Result<()> {
    writer.write_all(&[value])
}

fn write_u16<W: Write>(writer: &mut W, value: u16) -> io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

fn write_u32<W: Write>(writer: &mut W, value: u32) -> io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

fn write_u64<W: Write>(writer: &mut W, value: u64) -> io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

fn write_bytes<W: Write>(writer: &mut W, bytes: &[u8]) -> io::Result<()> {
    let len = u32::try_from(bytes.len()).map_err(|_| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("field length {} exceeded u32", bytes.len()),
        )
    })?;
    write_u32(writer, len)?;
    writer.write_all(bytes)
}

fn read_u8<R: Read>(reader: &mut R) -> io::Result<u8> {
    let mut buf = [0u8; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn read_u8_opt<R: Read>(reader: &mut R) -> io::Result<Option<u8>> {
    let mut buf = [0u8; 1];
    match reader.read_exact(&mut buf) {
        Ok(()) => Ok(Some(buf[0])),
        Err(err) if err.kind() == ErrorKind::UnexpectedEof => Ok(None),
        Err(err) => Err(err),
    }
}

fn read_u16<R: Read>(reader: &mut R) -> io::Result<u16> {
    let mut buf = [0u8; 2];
    reader.read_exact(&mut buf)?;
    Ok(u16::from_le_bytes(buf))
}

fn read_u32<R: Read>(reader: &mut R) -> io::Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_u32_opt<R: Read>(reader: &mut R) -> io::Result<Option<u32>> {
    let mut buf = [0u8; 4];
    match reader.read_exact(&mut buf) {
        Ok(()) => Ok(Some(u32::from_le_bytes(buf))),
        Err(err) if err.kind() == ErrorKind::UnexpectedEof => Ok(None),
        Err(err) => Err(err),
    }
}

fn read_u64<R: Read>(reader: &mut R) -> io::Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_bytes<R: Read>(reader: &mut R) -> io::Result<Vec<u8>> {
    let len = read_u32(reader)? as usize;
    let mut bytes = vec![0u8; len];
    reader.read_exact(&mut bytes)?;
    Ok(bytes)
}

fn read_strand<R: Read>(reader: &mut R) -> io::Result<ISOMSTRAND> {
    ISOMSTRAND::try_from(read_u8(reader)?)
        .map_err(|err| io::Error::new(ErrorKind::InvalidData, err.to_string()))
}

fn bytes_to_string(bytes: Vec<u8>) -> Result<String, GTFError> {
    String::from_utf8(bytes)
        .map_err(|err| GTFError::Io(Error::new(ErrorKind::InvalidData, err.to_string())))
}

fn bytes_to_string_io(bytes: Vec<u8>) -> io::Result<String> {
    String::from_utf8(bytes).map_err(|err| io::Error::new(ErrorKind::InvalidData, err.to_string()))
}

#[derive(Error, Debug)]
pub enum GTFError {
    #[error("Invalid GTF format")]
    InvalidGTFFormat { line_no: usize },

    #[error("GTF must contain at least one transcript record")]
    MissingTranscriptRecord,

    #[error(
        "GTF transcript/exon chromosome mismatch. Transcript-only seqids: {transcript_only:?}; exon-only seqids: {exon_only:?}"
    )]
    TranscriptExonChromMismatch {
        transcript_only: Vec<String>,
        exon_only: Vec<String>,
    },

    #[error(transparent)]
    Io(#[from] io::Error),
}
