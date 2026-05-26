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

pub struct RefPTIRManager {
    pub filename: PathBuf,
    pub ptirs: Vec<RefPTIR>,
    pub string_pool: StringPool,
    intervals_map: HashMap<String, Lapper<u32, usize>>,
}

impl RefPTIRManager {
    pub fn open<P: AsRef<Path>>(gtf_path: P) -> Result<Self, ClassifyError> {
        let gtf_path = gtf_path.as_ref().to_path_buf();
        let mut isomx_path = gtf_path.clone();
        isomx_path.add_extension("isomx");

        let f = File::open(isomx_path)?;
        let mut index_reader = IndexReader::open(f, 0)?;

        let mut attr_string_pool = StringPool::new();

        let attr_lines = extract_attr_lines(&gtf_path)?;

        let mut ptirs: Vec<RefPTIR> = Vec::new();
        let mut temp_intervals: HashMap<String, Vec<Interval<u32, usize>>> = HashMap::default();

        let mut chr_maps = index_reader.get_chromosome_readers_map()?;
        for (chr_name, chrom_block_builder) in &mut chr_maps {
            let chr_ivs = temp_intervals.entry(chr_name.clone()).or_default();
            while let Some(txbase) = chrom_block_builder.next_record()? {
                let ptir = PTIR::from_tx_base(
                    txbase,
                    0,
                    &chrom_block_builder.junction_pool,
                    &chrom_block_builder.splice_site_pool,
                    &chrom_block_builder.string_pool,
                );

                let attr_kvs = attr_lines.get(&ptir.source_txid).ok_or_else(|| {
                    ClassifyError::FailedParseGTF {
                        reason: format!(
                            "Transcript ID {} not found in reference GTF attributes.",
                            ptir.source_txid
                        ),
                    }
                })?;

                let gene_name = attr_kvs
                    .iter()
                    .find(|(k, _)| k == "gene_name")
                    .map(|(_, v)| v.clone())
                    .unwrap_or_default();

                let mut span_vec = Vec::new();
                for (k, v) in attr_kvs {
                    let k_span = attr_string_pool.add(k)?;
                    let v_span = attr_string_pool.add(v)?;
                    span_vec.push((k_span, v_span));
                }

                let ptir_idx = ptirs.len();
                chr_ivs.push(Interval {
                    start: ptir.start,
                    stop: ptir.end,
                    val: ptir_idx,
                });
                ptirs.push(RefPTIR {
                    base: ptir,
                    attrs: span_vec,
                    gene_name,
                });
            }
        }

        attr_string_pool.shrink_to_read_only();

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

    pub fn find_ovlp(&self, chr_name: &str, start: u32, end: u32) -> Option<Vec<&RefPTIR>> {
        let lapper = self.intervals_map.get(chr_name)?;
        let results: Vec<&RefPTIR> = lapper
            .find(start, end)
            .map(|iv| &self.ptirs[iv.val])
            .collect();
        if results.is_empty() {
            None
        } else {
            Some(results)
        }
    }
}

impl LogMemSize for RefPTIRManager {
    fn get_mem_size(&self) -> usize {
        use std::mem::size_of;

        let mut total = self.ptirs.capacity() * size_of::<RefPTIR>();

        for rp in &self.ptirs {
            total += rp.base.source_txid.capacity();
            total += rp.base.source_geneid.capacity();
            total += rp.gene_name.capacity();
            if let Some(jv) = &rp.base.junction_vec {
                total += jv.capacity() * size_of::<(u32, u32)>();
            }
            if let Some(sv) = &rp.base.splice_site_vec {
                total += sv.capacity() * size_of::<SpliceSitePair>();
            }
            total += rp.attrs.capacity() * size_of::<(StringSpan, StringSpan)>();
        }

        total += self.string_pool.heap_bytes();

        total +=
            self.intervals_map.capacity() * (size_of::<String>() + size_of::<Lapper<u32, usize>>());
        for (chr_name, lapper) in &self.intervals_map {
            total += chr_name.capacity();
            // intervals Vec
            total += lapper.intervals.capacity() * size_of::<Interval<u32, usize>>();
            // starts and stops Vecs are private but same length as intervals
            total += lapper.len() * 2 * size_of::<u32>();
        }

        total += self.filename.as_os_str().len();

        total
    }
}

fn breakdown_attrs(attr: &str) -> Vec<(String, String)> {
    let mut attrs = Vec::new();
    for part in attr.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let mut it = part.splitn(2, ' ');
        let k = match it.next() {
            Some(k) => k.trim(),
            None => continue,
        };
        let v = match it.next() {
            Some(v) => v.trim().trim_matches('"'),
            None => continue,
        };
        if !k.is_empty() {
            attrs.push((k.to_string(), v.to_string()));
        }
    }
    attrs
}

pub fn extract_attr_lines<P: AsRef<Path>>(
    p: P,
) -> Result<HashMap<String, Vec<(String, String)>>, ClassifyError> {
    let mut reader = crate::utils::open_file_bufread(p)?;
    let mut line = String::new();
    let mut out: HashMap<String, Vec<(String, String)>> = HashMap::default();

    while let Ok(bytes) = reader.read_line(&mut line) {
        if bytes == 0 {
            break;
        }
        let parts: Vec<&str> = line.splitn(9, '\t').collect();
        if parts.len() >= 9 && parts[2] == "transcript" {
            let transcript_id = parts[8]
                .split(';')
                .find_map(|kv| {
                    let mut it = kv.trim().splitn(2, ' ');
                    let k = it.next()?;
                    if k == "transcript_id" {
                        Some(it.next()?.trim().trim_matches('"').to_string())
                    } else {
                        None
                    }
                })
                .ok_or_else(|| ClassifyError::FailedParseGTF {
                    reason: "Can not read transcript_id from the reference GTF.".to_string(),
                })?;
            if out.contains_key(&transcript_id) {
                return Err(ClassifyError::FailedParseGTF {
                    reason: format!(
                        "Duplicate transcript_id {} found in reference GTF.",
                        transcript_id
                    ),
                });
            }
            out.insert(transcript_id, breakdown_attrs(parts[8]));
        }
        line.clear();
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_GTF: &str = "test/gencode.v49.basic.annotation.sorted.gtf.gz";

    #[test]
    fn open_loads_all_transcripts_and_attrs() {
        let mgr = RefPTIRManager::open(TEST_GTF).expect("RefPTIRManager::open failed");

        // total transcript count matches the GTF
        assert_eq!(mgr.ptirs.len(), 280000);

        // locate a known transcript
        let tx = mgr
            .ptirs
            .iter()
            .find(|p| p.base.source_txid == "ENST00000832828.1")
            .expect("ENST00000832828.1 not found");

        assert_eq!(tx.base.start, 11426);
        assert_eq!(tx.base.end, 14409);

        // gene_id attribute is present and correct
        let gene_id = tx
            .attrs
            .iter()
            .find_map(|(k_span, v_span)| {
                if mgr.string_pool.get(*k_span).ok()? == "gene_id" {
                    mgr.string_pool.get(*v_span).ok()
                } else {
                    None
                }
            })
            .expect("gene_id attr not found");
        assert_eq!(gene_id, "ENSG00000290825.2");
    }

    #[test]
    fn find_ovlp_returns_correct_hits() {
        let mgr = RefPTIRManager::open(TEST_GTF).expect("RefPTIRManager::open failed");

        // chr1:11000-12500 overlaps ENST00000832828.1 (11426-14409) and ENST00000450305.2 (12010-13670)
        let hits = mgr
            .find_ovlp("chr1", 11000, 12500)
            .expect("expected overlapping transcripts");
        assert_eq!(hits.len(), 2);
        let mut ids: Vec<&str> = hits.iter().map(|p| p.base.source_txid.as_str()).collect();
        ids.sort_unstable();
        assert_eq!(ids, ["ENST00000450305.2", "ENST00000832828.1"]);

        // empty region before all chr1 transcripts
        assert!(mgr.find_ovlp("chr1", 1, 100).is_none());

        // non-existent chromosome
        assert!(mgr.find_ovlp("chrZZZ", 0, 1_000_000).is_none());
    }

    #[test]
    fn get_mem_size_is_in_expected_range() {
        use crate::traits::LogMemSize;

        let mgr = RefPTIRManager::open(TEST_GTF).expect("RefPTIRManager::open failed");
        let bytes = mgr.get_mem_size();
        let mb = bytes as f64 / 1024.0 / 1024.0;
        eprintln!("RefPTIRManager::get_mem_size = {bytes} bytes ({mb:.1} MB)");

        // Vec<RefPTIR> grows by doubling: 280K elements land in a 524K-slot allocation
        // (524K × 208 B ≈ 109 MB) rather than the logical 280K × 208 B ≈ 58 MB.
        // Measured: ~240 MB on gencode.v49 (280K transcripts).
        assert!(
            bytes > 150_000_000,
            "too small — likely missing a heap term"
        );
        assert!(bytes < 350_000_000, "too large — likely double-counting");
    }
}
