use std::io::{Error, ErrorKind};
use std::{collections::VecDeque, io::BufRead, path::Path};

// use clap::error;
use log::error;
use log::warn;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::{core::tx_strand::ISOMSTRAND, utils::open_file_bufread};
use thiserror::Error;

/// Scan the GTF once to collect sorted unique chrom names and compute a
/// content hash (xxh3-128 of the decompressed bytes) and file size.
pub fn profile_gtf<P: AsRef<Path>>(path: P) -> Result<(Vec<String>, [u8; 16], u64), GTFError> {
    let file_size = std::fs::metadata(path.as_ref())?.len();

    let mut bufreader = open_file_bufread(path)?;
    let mut hasher = xxhash_rust::xxh3::Xxh3::new();
    let mut line = String::new();
    let mut transcript_chroms: FxHashSet<String> = FxHashSet::default();
    let mut exon_chroms: FxHashSet<String> = FxHashSet::default();
    let mut has_transcript = false;
    let mut line_no = 0usize;

    loop {
        let n = bufreader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        line_no += 1;
        hasher.update(line.as_bytes());

        if line.starts_with('#') {
            line.clear();
            continue;
        }

        let mut cols = line.split('\t');
        let chrom_name = cols
            .nth(0)
            .ok_or(GTFError::InvalidGTFFormat { line_no })?
            .to_string();
        cols.next(); // skip source
        let feature = cols.next().ok_or(GTFError::InvalidGTFFormat { line_no })?;

        if feature != "transcript" && feature != "exon" {
            line.clear();
            continue;
        }

        cols.next()
            .ok_or(GTFError::InvalidGTFFormat { line_no })?
            .parse::<u32>()
            .map_err(|_| GTFError::InvalidGTFFormat { line_no })?;

        if feature == "transcript" {
            has_transcript = true;
            transcript_chroms.insert(chrom_name);
        } else {
            exon_chroms.insert(chrom_name);
        }
        line.clear();
    }

    if !has_transcript {
        return Err(GTFError::MissingTranscriptRecord);
    }

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

    let mut chrom_names: Vec<String> = transcript_chroms.into_iter().collect();
    chrom_names.sort();

    let hash = hasher.digest128().to_le_bytes();
    Ok((chrom_names, hash, file_size))
}

pub enum GTFRecord {
    TxAttrs(TxAttrs),
    TxStructure(TxStructure),
}

#[derive(Debug, Clone)]
pub struct TxAttrs {
    chrname: String,
    attr_string: String,
}

impl TxAttrs {
    pub fn attr_string(&self) -> &str {
        &self.attr_string
    }

    pub fn chrname(&self) -> &str {
        &self.chrname
    }
}

/// GTF tx record
#[derive(Debug, Clone)]
pub struct TxStructure {
    pub gidx: u32,
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub strand: ISOMSTRAND,
    pub exons: Vec<(u32, u32)>,
    pub tx_id: String,
    pub gene_id: String,
    pub is_empty: bool,
}

impl TxStructure {
    pub fn default() -> Self {
        Self {
            gidx: 0,
            chrom: "".to_string(),
            start: 0,
            end: 0,
            strand: ISOMSTRAND::Unknown,
            exons: Vec::new(),
            tx_id: "".to_string(),
            gene_id: "".to_string(),
            is_empty: true,
        }
    }

    pub fn set_gidx(&mut self, idx: u32) {
        self.gidx = idx;
        // self.is_empty = false;
    }

    pub fn set_chrom(&mut self, chrom: String) {
        self.chrom = chrom;
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

    /// return the offset for all exons, base point starts from left position of the left most exon
    /// eg.  (10,20) , (30,40) -> (0,10), (20,30)
    pub fn get_0based_exon_relative_offset(&self) -> Vec<(u32, u32)> {
        let offsets = self.exons.clone();
        let base = offsets[0].0;
        let offsets: Vec<(u32, u32)> = offsets
            .into_iter()
            .map(|item| (item.0 - base, item.1 - base + 1))
            .collect();
        offsets
    }
}

pub struct MyGTFReader {
    pub ready_txs: VecDeque<TxStructure>,
}

impl MyGTFReader {
    pub fn new<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let mut bufreader = open_file_bufread(path)?;
        let mut txs = FxHashMap::default();
        let mut line = String::new();

        loop {
            line.clear();
            if bufreader.read_line(&mut line)? == 0 {
                break;
            }

            if line.starts_with('#') {
                continue;
            }

            let (chrom, feat, start, end, strand, tx_id, gene_id) = process_gtf_line(&line)?;
            if feat != "exon" {
                continue;
            }

            // validate the start and end coordinates of exons.
            // in case of invalide record has same start and end
            if start > end {
                warn!(
                    "Invalid GTF record with start > end, affected line: {}",
                    line
                );
                continue;
            }

            Self::add_exon_to_tx_map(&mut txs, chrom, start, end, strand, tx_id, gene_id);
        }

        let mut ready_txs: Vec<TxStructure> = txs
            .drain()
            .map(|(_, mut tx)| {
                tx.sort_exons();
                tx
            })
            .collect();

        ready_txs.sort_by(|a, b| {
            (
                a.chrom.as_str(),
                a.start,
                a.end,
                a.strand.clone(),
                a.tx_id.as_str(),
                a.gene_id.as_str(),
            )
                .cmp(&(
                    b.chrom.as_str(),
                    b.start,
                    b.end,
                    b.strand.clone(),
                    b.tx_id.as_str(),
                    b.gene_id.as_str(),
                ))
        });

        for (idx, tx) in ready_txs.iter_mut().enumerate() {
            tx.set_gidx(idx as u32);
        }

        Ok(Self {
            ready_txs: ready_txs.into(),
        })
    }

    fn add_exon_to_tx_map(
        txs: &mut FxHashMap<String, TxStructure>,
        chrom: String,
        start: u32,
        end: u32,
        strand: ISOMSTRAND,
        tx_id: String,
        gene_id: String,
    ) {
        let tx = txs.entry(tx_id.clone()).or_insert_with(|| {
            let mut tx = TxStructure::default();
            tx.set_start(start);
            tx.set_end(end);
            tx.set_chrom(chrom.clone());
            tx.set_strand(strand);
            tx.set_tx_id(tx_id.clone());
            tx.set_gene_id(gene_id.clone());
            tx
        });

        if tx.gene_id != gene_id {
            warn!(
                "Transcript {} has inconsistent gene_idin it's exon record. {}: {} vs {}",
                tx_id, chrom, tx.gene_id, gene_id
            );
        }

        if tx.strand != strand {
            warn!(
                "Transcript {} has inconsistent strand in it's exon record at chr: {}. {} vs {}",
                tx_id, chrom, tx.strand, strand
            );
        }

        tx.add_exon((start, end));
    }

    pub fn next(&mut self) -> Result<Option<TxStructure>, GTFError> {
        Ok(self.ready_txs.pop_front())
    }
}

/// process one line of GTF file, return chrom, feature type, start, end, strand,
/// transcript_id and gene_id. The start and end are 1-based and end is inclusive.
pub fn process_gtf_line(
    s: &str,
) -> Result<
    (
        String,     // chrom (col 0)
        String,     // feature_type (col 2): "transcript" / "exon" / ...
        u32,        // start (1-based)
        u32,        // end   (1-based, inclusive)
        ISOMSTRAND, // strand: 0=+, 1=-
        String,     // transcript_id
        String,     // gene_id
    ),
    Error,
> {
    let parts: Vec<&str> = s.split('\t').collect();

    if parts.len() < 9 {
        return Err(Error::new(
            ErrorKind::InvalidData,
            format!("Invalid GTF line: fewer than 9 columns. Affected line: {}", s.trim_end()),
        ));
    }

    // col 0: chrom
    let chrom = parts[0].to_string();
    // col 2: feature type
    let feature_type = parts[2].to_string();
    // col 3: start, 1-based
    let start = parts[3].parse::<u32>().expect("invalid start coordinate");
    // col 4: end, 1-based inclusive
    let end = parts[4].parse::<u32>().expect("invalid end coordinate");
    // col 6: strand
    let strand = match parts[6] {
        "-" => ISOMSTRAND::Minus,
        "+" => ISOMSTRAND::Plus,
        "." => ISOMSTRAND::Unknown,
        _ => {
            error!("Unknown Strand for transcript: {s} ");
            std::process::exit(1);
        }
    };
    // col 8: attributes
    let (tx_id, gene_id) = parse_gtf_attributes(parts[8]);

    if tx_id.is_empty() || gene_id.is_empty() {
        let missing = match (tx_id.is_empty(), gene_id.is_empty()) {
            (true, true) => "transcript_id and gene_id",
            (true, false) => "transcript_id",
            (false, true) => "gene_id",
            (false, false) => unreachable!(),
        };
        return Err(Error::new(
            ErrorKind::InvalidData,
            format!(
                "Missing required GTF attribute(s): {missing}. Affected line: {}",
                s.trim_end()
            ),
        ));
    }

    Ok((chrom, feature_type, start, end, strand, tx_id, gene_id))
}

/// take the attributes column of a GTF line and extract the transcript_id and gene_id values, supporting both quoted and unquoted formats.
fn parse_gtf_attributes(attrs: &str) -> (String, String) {
    let mut tx_id = String::new();
    let mut gene_id = String::new();

    for attr in attrs.split(';') {
        let attr = attr.trim();
        if attr.is_empty() {
            continue;
        }
        if attr.starts_with("transcript_id") {
            tx_id = extract_attr_value(attr);
        } else if attr.starts_with("gene_id") {
            gene_id = extract_attr_value(attr);
        }
        if !tx_id.is_empty() && !gene_id.is_empty() {
            break; // 两个都找到了，提前退出
        }
    }

    (tx_id, gene_id)
}

/// extract the value part of a GTF attribute, supporting both quoted and unquoted formats.
fn extract_attr_value(attr: &str) -> String {
    // quoted: gene_id "ENSG00000000003";
    if let Some(q_start) = attr.find('"') {
        if let Some(q_len) = attr[q_start + 1..].find('"') {
            return attr[q_start + 1..q_start + 1 + q_len].to_string();
        }
    }
    // unquoted: gene_id ENSG00000000003
    attr.split_ascii_whitespace()
        .nth(1)
        .unwrap_or("")
        .to_string()
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
    Io(#[from] std::io::Error),
}
