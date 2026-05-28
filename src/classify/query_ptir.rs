//! struct of query PIIR

use std::{
    fs::File,
    path::{Path, PathBuf},
};

use crate::{
    classify::classify_error::ClassifyError,
    core::{ptir::PTIR, tx_strand::ISOMSTRAND, tx_type::TxType},
    index::{
        attributes_index::AttrIndexReader,
        reader::{ChromBlockReader, IndexReader},
    },
};

#[derive(Debug, Clone)]
pub struct QueryPTIR {
    pub chr_name: String,
    pub base: PTIR,
    pub attr_raw_string: Vec<u8>,
}

impl QueryPTIR {
    pub fn new(chr_name: &str, ptir: PTIR, attr_string: Vec<u8>) -> Self {
        Self {
            chr_name: chr_name.to_string(),
            base: ptir,
            attr_raw_string: attr_string,
        }
    }

    pub fn start(&self) -> u32 {
        self.base.start
    }

    pub fn end(&self) -> u32 {
        self.base.end
    }

    pub fn standard(&self) -> &ISOMSTRAND {
        &self.base.strand
    }

    pub fn n_exons(&self) -> u16 {
        self.base.n_exons
    }

    pub fn junction_vec(&self) -> &Option<Vec<(u32, u32)>> {
        &self.base.junction_vec
    }

    pub fn tx_type(&self) -> &TxType {
        &self.base.tx_type
    }

    pub fn exons_vec(&self) -> Vec<(u32, u32)> {
        let Some(junctions) = self.base.junction_vec.as_ref() else {
            return vec![(self.start(), self.end())];
        };

        if junctions.is_empty() {
            return vec![(self.start(), self.end())];
        }

        let mut boundaries = Vec::with_capacity(junctions.len() * 2 + 2);
        boundaries.push(self.start());
        boundaries.extend(junctions.iter().flat_map(|&(left, right)| [left, right]));
        boundaries.push(self.end());

        boundaries
            .chunks_exact(2)
            .map(|chunk| (chunk[0], chunk[1]))
            .collect()
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
        let mut attr_file_name = gtf_path.clone();
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

    pub fn next_record(&mut self) -> Option<QueryPTIR> {
        loop {
            let txbase = match self.current_reader.next_record() {
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
            return Some(QueryPTIR::new(
                &self.current_reader.chrom_name,
                ptir,
                attr_bytes,
            ));
        }
    }
}
