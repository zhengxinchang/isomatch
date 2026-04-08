use std::{collections::VecDeque, io::BufRead, path::Path};

use log::warn;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::utils::open_file_bufread;
use thiserror::Error;

/// Scan the GTF once to collect ordered unique chrom names (transcript lines only).
/// Returns an error if the GTF is not sorted by chrom then by start position.
pub fn profile_gtf<P: AsRef<std::path::Path>>(path: P) -> Result<Vec<String>, GTFError> {
    let mut bufreader = open_file_bufread(path).map_err(|_| GTFError::IoError {})?;
    let mut line = String::new();
    let mut chrom_set: FxHashSet<String> = FxHashSet::default();
    let mut chrom_names: Vec<String> = Vec::new();
    let mut prev_chrom = String::new();
    let mut prev_start = 0u32;
    let mut line_no = 0usize;

    while let Ok(n) = bufreader.read_line(&mut line) {
        if n == 0 {
            break;
        }
        line_no += 1;

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

        if feature != "transcript" {
            line.clear();
            continue;
        }

        let start = cols
            .next()
            .ok_or(GTFError::InvalidGTFFormat { line_no })?
            .parse::<u32>()
            .map_err(|_| GTFError::InvalidGTFFormat { line_no })?;

        if !prev_chrom.is_empty() && chrom_name != prev_chrom {
            if chrom_set.contains(&chrom_name) {
                return Err(GTFError::UnsortedGTF { line: line });
            }
            prev_start = 0;
        }
        if chrom_name == prev_chrom && prev_start > start {
            return Err(GTFError::UnsortedGTF { line: line });
        }

        if !chrom_set.contains(&chrom_name) {
            chrom_set.insert(chrom_name.clone());
            chrom_names.push(chrom_name.clone());
        }
        prev_chrom = chrom_name;
        prev_start = start;
        line.clear();
    }

    Ok(chrom_names)
}

/// GTF tx record, both for input GTF and referencce annotation GTF.
#[derive(Debug, Clone)]
pub struct GTFTx {
    pub idx: u32,
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub strand: u8,
    pub exons: Vec<(u32, u32)>,
    pub tx_id: String,
    pub gene_id: String,
    pub is_empty: bool,
}

impl GTFTx {
    pub fn default() -> Self {
        Self {
            idx: 0,
            chrom: "".to_string(),
            start: 0,
            end: 0,
            strand: 0,
            exons: Vec::new(),
            tx_id: "".to_string(),
            gene_id: "".to_string(),
            is_empty: true,
        }
    }

    pub fn set_idx(&mut self, idx: u32) {
        self.idx = idx;
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

    pub fn set_strand(&mut self, strand: u8) {
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
    pub bufreader: Box<dyn BufRead>,
    pub current_tx_idx: u32,
    pub current_chrom: Option<String>,
    pub chrom_txs: FxHashMap<String, GTFTx>,
    pub ready_txs: VecDeque<GTFTx>,
}

impl MyGTFReader {
    pub fn new<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        Ok(Self {
            bufreader: open_file_bufread(path)?,
            current_tx_idx: 0,
            current_chrom: None,
            chrom_txs: FxHashMap::default(),
            ready_txs: VecDeque::new(),
        })
    }

    fn add_exon_to_current_chrom(
        &mut self,
        chrom: String,
        start: u32,
        end: u32,
        strand: u8,
        tx_id: String,
        gene_id: String,
    ) {
        let tx = self.chrom_txs.entry(tx_id.clone()).or_insert_with(|| {
            let mut tx = GTFTx::default();
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
                "Transcript {} has inconsistent gene_id within chromosome {}: {} vs {}",
                tx_id, chrom, tx.gene_id, gene_id
            );
        }

        if tx.strand != strand {
            warn!(
                "Transcript {} has inconsistent strand within chromosome {}: {} vs {}",
                tx_id, chrom, tx.strand, strand
            );
        }

        tx.add_exon((start, end));
    }

    fn flush_current_chrom(&mut self) {
        if self.chrom_txs.is_empty() {
            self.current_chrom = None;
            return;
        }

        let mut txs: Vec<GTFTx> = self
            .chrom_txs
            .drain()
            .map(|(_, mut tx)| {
                tx.sort_exons();
                tx
            })
            .collect();

        txs.sort_by(|a, b| {
            (
                a.start,
                a.end,
                a.strand,
                a.tx_id.as_str(),
                a.gene_id.as_str(),
            )
                .cmp(&(
                    b.start,
                    b.end,
                    b.strand,
                    b.tx_id.as_str(),
                    b.gene_id.as_str(),
                ))
        });

        for mut tx in txs {
            tx.set_idx(self.current_tx_idx);
            self.current_tx_idx += 1;
            self.ready_txs.push_back(tx);
        }

        self.current_chrom = None;
    }
}

impl Iterator for MyGTFReader {
    type Item = GTFTx;

    /// if the ready_txs is empty, this
    /// next() will load all the records in next chromsome,
    /// and rebuild the ready_txs.
    /// this require the gtf sorted by at least the chromosome.
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(tx) = self.ready_txs.pop_front() {
            return Some(tx);
        }

        let mut line = String::new();
        loop {
            line.clear();
            let n = match self.bufreader.read_line(&mut line) {
                Ok(n) => n,
                Err(_) => {
                    self.flush_current_chrom();
                    return self.ready_txs.pop_front();
                }
            };

            if n == 0 {
                self.flush_current_chrom();
                return self.ready_txs.pop_front();
            }

            if line.starts_with('#') {
                continue;
            }

            let (chrom, feat, start, end, strand, tx_id, gene_id) = process_gtf_line(&line);

            // only process exon lines
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

            if let Some(current_chrom) = &self.current_chrom {
                if current_chrom != &chrom {
                    self.flush_current_chrom();
                }
            }

            if self.current_chrom.is_none() {
                self.current_chrom = Some(chrom.clone());
            }

            self.add_exon_to_current_chrom(chrom, start, end, strand, tx_id, gene_id);

            if let Some(tx) = self.ready_txs.pop_front() {
                return Some(tx);
            }
        }
    }
}

/// process one line of GTF file, return chrom, feature type, start, end, strand,
/// transcript_id and gene_id. The start and end are 1-based and end is inclusive.
pub fn process_gtf_line(
    s: &str,
) -> (
    String, // chrom (col 0)
    String, // feature_type (col 2): "transcript" / "exon" / ...
    u32,    // start (1-based)
    u32,    // end   (1-based, inclusive)
    u8,     // strand: 0=+, 1=-
    String, // transcript_id
    String, // gene_id
) {
    let parts: Vec<&str> = s.split('\t').collect();

    if parts.len() < 9 {
        panic!(
            "Invalid GTF line: fewer than 9 columns, affected line: {}",
            s
        );
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
        "-" => 1u8,
        _ => 0u8, // '+' or '.' → forward
    };
    // col 8: attributes
    let (tx_id, gene_id) = parse_gtf_attributes(parts[8]);

    (chrom, feature_type, start, end, strand, tx_id, gene_id)
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
    #[error("Can not read GTF file")]
    IoError {},

    #[error("Unsorted GTF: {line}")]
    UnsortedGTF { line: String },

    #[error("Invalid GTF format")]
    InvalidGTFFormat { line_no: usize },
}

#[cfg(test)]
mod tests {
    use super::MyGTFReader;
    use flate2::{Compression, write::GzEncoder};
    use std::{
        fs::{self, File},
        io::Write,
        path::{Path, PathBuf},
        time::{SystemTime, UNIX_EPOCH},
    };

    fn unique_temp_path(ext: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "isomatch_gtf_{}_{}.{}",
            std::process::id(),
            nanos,
            ext
        ))
    }

    fn write_plain_gtf(path: &Path, content: &str) {
        fs::write(path, content).unwrap();
    }

    fn write_gzip_gtf(path: &Path, content: &str) {
        let file = File::create(path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(content.as_bytes()).unwrap();
        encoder.finish().unwrap();
    }

    #[test]
    fn my_gtf_reader_supports_plain_text_gtf() {
        let path = unique_temp_path("gtf");
        let content = concat!(
            "chr1\tsrc\texon\t1\t10\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n",
            "chr1\tsrc\texon\t21\t30\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n",
            "chr1\tsrc\texon\t41\t50\t.\t+\t.\tgene_id \"g2\"; transcript_id \"tx2\";\n",
        );
        write_plain_gtf(&path, content);

        let records: Vec<_> = MyGTFReader::new(&path).unwrap().collect();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].tx_id, "tx1");
        assert_eq!(records[0].exons, vec![(1, 10), (21, 30)]);
        assert_eq!(records[1].tx_id, "tx2");

        fs::remove_file(path).unwrap();
    }

    #[test]
    fn my_gtf_reader_supports_gzip_gtf() {
        let path = unique_temp_path("gtf.gz");
        let content = concat!(
            "chr1\tsrc\texon\t1\t10\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n",
            "chr1\tsrc\texon\t21\t30\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n",
            "chr1\tsrc\texon\t41\t50\t.\t+\t.\tgene_id \"g2\"; transcript_id \"tx2\";\n",
        );
        write_gzip_gtf(&path, content);

        let records: Vec<_> = MyGTFReader::new(&path).unwrap().collect();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].tx_id, "tx1");
        assert_eq!(records[0].exons, vec![(1, 10), (21, 30)]);
        assert_eq!(records[1].tx_id, "tx2");

        fs::remove_file(path).unwrap();
    }

    #[test]
    fn my_gtf_reader_merges_interleaved_exons_within_chromosome() {
        let path = unique_temp_path("gtf");
        let content = concat!(
            "chr1\tsrc\texon\t1\t10\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n",
            "chr1\tsrc\texon\t11\t20\t.\t+\t.\tgene_id \"g2\"; transcript_id \"tx2\";\n",
            "chr1\tsrc\texon\t41\t50\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n",
            "chr1\tsrc\texon\t31\t40\t.\t+\t.\tgene_id \"g2\"; transcript_id \"tx2\";\n",
        );
        write_plain_gtf(&path, content);

        let records: Vec<_> = MyGTFReader::new(&path).unwrap().collect();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].tx_id, "tx1");
        assert_eq!(records[0].exons, vec![(1, 10), (41, 50)]);
        assert_eq!(records[1].tx_id, "tx2");
        assert_eq!(records[1].exons, vec![(11, 20), (31, 40)]);

        fs::remove_file(path).unwrap();
    }
}
