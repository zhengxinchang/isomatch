use std::{
    fs::File,
    io::BufRead,
    path::{Path, PathBuf},
};

use ahash::HashMap;
use rust_lapper::{Interval, Lapper};

use crate::{
    classify::classify_error::ClassifyError,
    core::{
        ptir::PTIR,
        splice_site_pair::SpliceSitePair,
        string_pool::{StringPool, StringSpan},
        tx_strand::ISOMSTRAND,
        tx_type::TxType,
    },
    index::reader::IndexReader,
    traits::LogMemSize,
};

pub struct RefPTIR {
    pub base: PTIR,
    pub attrs: Vec<(StringSpan, StringSpan)>, // key,value
    pub gene_name: String,
}

impl RefPTIR {
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
}
