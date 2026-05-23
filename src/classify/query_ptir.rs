//! struct of query PIIR

use std::path::{Path, PathBuf};

use crate::{
    core::ptir::PTIR,
    index::{
        attributes_index::{AttrIndexReader, RawStringSpan},
        reader::IndexReader,
    },
};

pub struct QueryPTIR {
    pub is_isomatch_merged: bool,
    pub base: PTIR,
    pub attr_raw_string: Vec<u8>,
}

impl QueryPTIR {
    pub fn new(ptir: &PTIR, attr_string: Vec<u8>) -> Self {
        // construct the QueryPTIR
        todo!()
    }
}

pub struct QueryPTIRManager {
    index_file_name: PathBuf,
    attr_file_name: PathBuf,
    index_reader: IndexReader,
    attr_index_reader: AttrIndexReader,
    total_tx_n: usize, // from IndexReader
}

impl QueryPTIRManager {
    pub fn open<P: AsRef<Path>>(gtf_path: P) -> Self {
        todo!()
    }

    pub fn next(&mut self) -> Option<QueryPTIR> {
        todo!()
    }
}
