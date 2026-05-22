use std::{
    fs::File,
    io::{BufRead, Read},
    path::{Path, PathBuf},
};

use ahash::HashMap;
use rust_lapper::{Interval, Lapper};

use crate::{
    classify::classify_error::ClassifyError,
    core::{
        ptir::PTIR,
        string_pool::{StringPool, StringSpan},
        tx_base::{TxBase, TxBaseTrait},
    },
    index::{format::ChromBlockBuilder, reader::IndexReader}, merge::grouped_ptirs,
};

/// PTIR representation of reference GTF

pub struct RefPTIR {
    txbase: PTIR,
    attrs: Vec<(StringSpan, StringSpan)>, // key,value
}

pub struct RefPTIRManager {
    filename: PathBuf,
    ptirs: Vec<RefPTIR>,
    string_pool: StringPool,
    intervals_map: HashMap<String, Lapper<u32, usize>>,
}

impl RefPTIRManager {
    pub fn open<P: AsRef<Path>>(gtf_path: P) -> Result<Self, ClassifyError> {
        let gtf_path = gtf_path.as_ref().to_path_buf();
        let mut isomx_path = gtf_path.clone();
        isomx_path.set_extension("isomx");

        let f = File::open(isomx_path)?;
        let mut index_reader = IndexReader::open(f, 0)?;

        let mut attr_string_pool = StringPool::new();
        let attr_lines = extract_attr_lines(&gtf_path)?;

        let mut ptirs: Vec<RefPTIR> = Vec::new();
        let mut temp_intervals: HashMap<String, Vec<Interval<u32, usize>>> = HashMap::default();

        let mut chr_maps = index_reader.get_chromosome_readers_map()?;
        for (chr_name, chrom_block_builder) in &mut chr_maps {
            while let Some(txbase) = chrom_block_builder.next()? {
                let ptir = PTIR::from_tx_base(
                    txbase,
                    0,
                    &chrom_block_builder.junction_pool,
                    &chrom_block_builder.splice_site_pool,
                    &chrom_block_builder.string_pool,
                );

                let attr_kvs = breakdown_attrs(
                    attr_lines
                        .get(&ptir.source_txid)
                        .ok_or_else(|| ClassifyError::FailedParseGTF {
                            reason: format!("Transcript ID {} not found in reference GTF attributes.", ptir.source_txid),
                        })?,
                )?;

                let mut span_vec = Vec::new();
                for (k, v) in attr_kvs {
                    let k_span = attr_string_pool.add(&k)?;
                    let v_span = attr_string_pool.add(&v)?;
                    span_vec.push((k_span, v_span));
                }

                let ptir_idx = ptirs.len();
                let iv = Interval { start: ptir.start, stop: ptir.end, val: ptir_idx };

                ptirs.push(RefPTIR { txbase: ptir, attrs: span_vec });

                temp_intervals
                    .entry(chr_name.clone())
                    .or_default()
                    .push(iv);
            }
        }

        let intervals_map = temp_intervals
            .into_iter()
            .map(|(chr, ivs)| (chr, Lapper::new(ivs)))
            .collect();

        Ok(Self {
            filename: gtf_path,
            ptirs,
            string_pool: attr_string_pool,
            intervals_map,
        })
    }

    pub fn find_ovlp(&mut self, chr_name: &str, start: u32, end: u32) -> Option<Vec<&RefPTIR>> {
        let lapper = self.intervals_map.get_mut(chr_name)?;
        let results: Vec<&RefPTIR> = lapper
            .find(start, end)
            .map(|iv| &self.ptirs[iv.val])
            .collect();
        if results.is_empty() { None } else { Some(results) }
    }
}

pub fn breakdown_attrs(attr: &str) -> Result<Vec<(String, String)>, ClassifyError> {
    let mut attrs = Vec::new();
    let parts: Vec<&str> = attr.split(';').collect();
    for part in parts {
        let part2 = part.trim();
        let kv: Vec<&str> = part2.splitn(2, ' ').collect();

        let v = kv[1].trim_matches('"');
        let k = kv[0].trim();

        attrs.push((k.to_string(), v.to_string()));
    }
    Ok(attrs)
}

pub fn extract_attr_lines<P: AsRef<Path>>(p: P) -> Result<HashMap<String,String>, ClassifyError> {
    let mut reader = crate::utils::open_file_bufread(p)?;
    let mut line = String::new();
    let mut out = HashMap::default();

    while let Ok(bytes) = reader.read_line(&mut line) {
        if bytes == 0 {
            break;
        }
        let parts: Vec<&str> = line.splitn(9, '\t').collect();
        if parts.len() >= 9 && parts[2] == "transcript" {
            // find the transcript_id in the attributes column (9th column)
            let transcript_id = parts[8]
                .split(';')
                .find_map(|kv| {
                    let kv2: Vec<&str> = kv.trim().splitn(2, ' ').collect();
                    if kv2.len() == 2 && kv2[0] == "transcript_id" {
                        Some(kv2[1].trim_matches('"').to_string())
                    } else {
                        None
                    }
                })
                .ok_or_else(|| ClassifyError::FailedParseGTF { reason: "Can not read transcript_id from the reference GTF.".to_string() } )?;
            if out.contains_key(&transcript_id){
                return Err(ClassifyError::FailedParseGTF { reason: format!("Duplicate transcript_id {} found in the reference GTF. Please ensure all transcript_id values are unique.", transcript_id) });
            }
            out.insert(transcript_id, parts[8].to_string());
        }
        line.clear();
    }
    Ok(out)
}
