//! struct of query PIIR

use std::{
    fs::File,
    path::{Path, PathBuf},
};

use crate::{
    classify::classify_error::ClassifyError,
    core::ptir::PTIR,
    index::{
        attributes_index::AttrIndexReader,
        reader::{ChromBlockReader, IndexReader},
    },
};

pub struct QueryPTIR {
    pub base: PTIR,
    pub attr_raw_string: Vec<u8>,
}

impl QueryPTIR {
    pub fn new(ptir: PTIR, attr_string: Vec<u8>) -> Self {
        Self {
            base: ptir,
            attr_raw_string: attr_string,
        }
    }
}

pub struct QueryPTIRManager {
    index_file_name: PathBuf,
    attr_file_name: PathBuf,
    index_reader: IndexReader,
    attr_index_reader: AttrIndexReader,
    total_tx_n: usize,
    index_chrnames: Vec<String>,
    current_reader: ChromBlockReader,
    chrom_idx: usize,
}

impl QueryPTIRManager {
    pub fn open<P: AsRef<Path>>(gtf_path: P) -> Result<Self, ClassifyError> {
        let gtf_path = gtf_path.as_ref().to_path_buf();
        let mut index_file_name = gtf_path.clone();
        index_file_name.add_extension("isomx");
        let mut attr_file_name = index_file_name.clone();
        attr_file_name.add_extension("isoms");

        let mut index_reader = IndexReader::open(File::open(&index_file_name)?, 0)?;
        let attr_index_reader = AttrIndexReader::open(&attr_file_name)?;

        let total_tx_n = index_reader.header.total_tx_n as usize;
        let index_chrnames = index_reader.chrom_names.clone();

        let first_chrom = index_chrnames.first().ok_or_else(|| {
            ClassifyError::Io(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "index contains no chromosomes",
            ))
        })?;
        let current_reader = index_reader.get_chromosome_reader(first_chrom)?;

        Ok(Self {
            index_file_name,
            attr_file_name,
            index_reader,
            attr_index_reader,
            total_tx_n,
            index_chrnames,
            current_reader,
            chrom_idx: 1,
        })
    }

    pub fn total_tx_n(&self) -> usize {
        self.total_tx_n
    }

    pub fn next(&mut self) -> Option<QueryPTIR> {
        loop {
            let txbase = match ChromBlockReader::next(&mut self.current_reader) {
                Ok(Some(tb)) => tb,
                Ok(None) => {
                    if self.chrom_idx >= self.index_chrnames.len() {
                        return None;
                    }
                    let next_chrom = self.index_chrnames[self.chrom_idx].clone();
                    self.current_reader = self
                        .index_reader
                        .get_chromosome_reader(&next_chrom)
                        .unwrap_or_else(|e| {
                            eprintln!("error loading chromosome {}: {}", next_chrom, e);
                            std::process::exit(1);
                        });
                    self.chrom_idx += 1;
                    continue;
                }
                Err(e) => {
                    eprintln!("error reading next transcript from query index: {}", e);
                    std::process::exit(1);
                }
            };
            let tx_gidx = txbase.tx_idx;
            let ptir = PTIR::from_tx_base(
                txbase,
                0,
                &self.current_reader.junction_pool,
                &self.current_reader.splice_site_pool,
                &self.current_reader.string_pool,
            );
            let attr_bytes = self
                .attr_index_reader
                .get_attr(tx_gidx)
                .unwrap_or(None)
                .unwrap_or_default();
            return Some(QueryPTIR::new(ptir, attr_bytes));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Instant;

    const TEST_GTF: &str = "test/gencode.v49.basic.annotation.sorted.gtf.gz";

    #[test]
    fn iterate_all_transcripts_timed() {
        let mut mgr = QueryPTIRManager::open(TEST_GTF).expect("QueryPTIRManager::open failed");

        let expected = mgr.total_tx_n();
        eprintln!("header.total_tx_n = {}", expected);

        let t0 = Instant::now();
        let mut count = 0usize;
        let mut first: Option<(String, u32, u32)> = None;
        while let Some(qtx) = mgr.next() {
            if first.is_none() {
                first = Some((qtx.base.source_txid.clone(), qtx.base.start, qtx.base.end));
            }
            count += 1;
        }
        let elapsed = t0.elapsed();

        eprintln!(
            "iterated {} transcripts in {:.3}s ({:.0} tx/s)",
            count,
            elapsed.as_secs_f64(),
            count as f64 / elapsed.as_secs_f64(),
        );
        if let Some((tx_id, start, end)) = first {
            eprintln!("first transcript: {} {}–{}", tx_id, start, end);
        }

        assert_eq!(count, 280000, "expected 280000 transcripts");
    }
}
