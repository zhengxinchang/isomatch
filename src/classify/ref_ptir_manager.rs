use std::{
    fs::File,
    path::{Path, PathBuf},
};

use ahash::{HashMap, HashSet};
use log::error;
use log::warn;
use rust_lapper::{Interval, Lapper};

use crate::{
    classify::{classify_error::ClassifyError, ref_ptir::RefPTIR},
    core::{
        ptir::PTIR,
        splice_site_pair::SpliceSitePair,
        string_pool::{StringPool, StringSpan},
        tx_strand::ISOMSTRAND,
    },
    index::reader::IndexReader,
    traits::LogMemSize,
};

type ChromId = u32;
type GeneId = u32;
type TxId = u32;
type Pos = u32;

#[derive(Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct Junction {
    start: Pos,
    end: Pos,
}

#[derive(Debug, Default)]
struct StringID {
    gene_to_id: HashMap<Box<str>, GeneId>,
    id_to_gene: Vec<Box<str>>,
    chrom_to_id: HashMap<Box<str>, ChromId>,
    id_to_chrom: Vec<Box<str>>,
}

impl StringID {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn get_or_insert_gene(&mut self, gene: &str) -> GeneId {
        if let Some(&id) = self.gene_to_id.get(gene) {
            return id;
        }

        let id = self.id_to_gene.len() as GeneId;
        let key: Box<str> = gene.into();

        self.id_to_gene.push(key.clone());
        self.gene_to_id.insert(key, id);

        id
    }

    pub fn get_or_insert_chrom(&mut self, chrom: &str) -> ChromId {
        if let Some(&id) = self.chrom_to_id.get(chrom) {
            return id;
        }

        let id = self.id_to_chrom.len() as ChromId;
        let key: Box<str> = chrom.into();

        self.id_to_chrom.push(key.clone());
        self.chrom_to_id.insert(key, id);

        id
    }

    pub fn gene_id(&self, gene: &str) -> Option<GeneId> {
        self.gene_to_id.get(gene).copied()
    }

    pub fn chrom_id(&self, chrom: &str) -> Option<ChromId> {
        self.chrom_to_id.get(chrom).copied()
    }

    pub fn gene_name(&self, id: GeneId) -> Option<&str> {
        self.id_to_gene.get(id as usize).map(|s| s.as_ref())
    }

    pub fn chrom_name(&self, id: ChromId) -> Option<&str> {
        self.id_to_chrom.get(id as usize).map(|s| s.as_ref())
    }

    pub fn gene_count(&self) -> usize {
        self.id_to_gene.len()
    }

    pub fn chrom_count(&self) -> usize {
        self.id_to_chrom.len()
    }
}

struct ChromIndex {
    starts: Vec<Pos>,
    ends: Vec<Pos>,

    // 分 strand 存 known junction pair
    junctions_plus: Vec<Junction>,
    junctions_minus: Vec<Junction>,

    // 用于 genic_intron: known intron intervals
    prefix_max_starts_plus: Vec<Pos>,
    prefix_max_starts_minus: Vec<Pos>,

    // 用于先找 overlap ref transcript
    refs_monoexon: Lapper<u32, usize>,
    refs_multiexon: Lapper<u32, usize>,
}
impl ChromIndex {
    pub fn new(
        mut junctions_plus: Vec<Junction>,
        mut junctions_minus: Vec<Junction>,
        interval_mono_exon: Vec<Interval<u32, usize>>,
        interval_multi_exon: Vec<Interval<u32, usize>>,
    ) -> Self {
        junctions_plus.sort();
        junctions_minus.sort();
        let mut starts = HashSet::default();
        let mut ends = HashSet::default();

        let mut prefix_max_starts_plus = Vec::new();
        let mut prefix_max_starts_minus = Vec::new();

        for junction in junctions_plus.iter() {
            starts.insert(junction.start);
            ends.insert(junction.end);
            prefix_max_starts_plus.push(junction.end)
        }

        for junction in junctions_minus.iter() {
            starts.insert(junction.start);
            ends.insert(junction.end);
            prefix_max_starts_minus.push(junction.end);
        }

        let mut starts: Vec<Pos> = starts.into_iter().collect();
        starts.sort();

        let mut ends: Vec<Pos> = ends.into_iter().collect();
        ends.sort();

        Self {
            starts,
            ends,
            junctions_plus,
            junctions_minus,
            prefix_max_starts_plus,
            prefix_max_starts_minus,
            refs_monoexon: Lapper::new(interval_mono_exon),
            refs_multiexon: Lapper::new(interval_multi_exon),
        }
    }
}

struct GeneIndex {
    junctions: Vec<Junction>,
    starts: Vec<Pos>,
    ends: Vec<Pos>,
}

impl GeneIndex {
    pub fn new(junctions: Option<&[(u32, u32)]>, start: u32, end: u32) -> Self {
        let mut juncs = Vec::new();
        if let Some(junctions) = junctions {
            for junc in junctions {
                juncs.push(Junction {
                    start: junc.0,
                    end: junc.1,
                })
            }
        }

        Self {
            junctions: juncs,
            starts: vec![start],
            ends: vec![end],
        }
    }
    pub fn add(&mut self, junctions: Option<&[(u32, u32)]>, start: u32, end: u32) {
        if let Some(junctions) = junctions {
            for junc in junctions {
                self.junctions.push(Junction {
                    start: junc.0,
                    end: junc.1,
                })
            }
        }

        self.starts.push(start);
        self.ends.push(start);
    }
}

pub struct RefPTIRManager {
    pub filename: PathBuf,
    pub ptirs: Vec<RefPTIR>,
    pub attr_string_pool: StringPool,
    // intervals_map: HashMap<String, Lapper<u32, usize>>,
    chroms: Vec<ChromIndex>,
    genes: Vec<GeneIndex>,
    stringids: StringID,
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
        // let mut temp_intervals: HashMap<String, Vec<Interval<u32, usize>>> = HashMap::default();

        let mut chr_maps = index_reader.get_chromosome_readers_map()?;

        let mut chrom_indexes = Vec::new();
        let mut geneid_indexes = Vec::new();
        let mut global_string_id = StringID::new();
        for (chr_name, chrom_block_builder) in &mut chr_maps {
            let mut interval_mono_exon: Vec<Interval<u32, usize>> = Vec::new();
            let mut interval_multi_exon: Vec<Interval<u32, usize>> = Vec::new();
            let mut junction_plus = Vec::new();
            let mut junction_minus = Vec::new();
            let chrom_id = global_string_id.get_or_insert_chrom(chr_name);
            while let Some(txbase) = chrom_block_builder.next_record()? {
                let ptir = PTIR::from_tx_base(
                    txbase,
                    0,
                    &chrom_block_builder.junction_pool,
                    &chrom_block_builder.splice_site_pool,
                    &chrom_block_builder.string_pool,
                );

                if matches!(ptir.strand, ISOMSTRAND::Unknown) {
                    warn!(
                        "Reference GTF contains unstranded transcript: {}",
                        &ptir.source_txid
                    );
                    continue;
                }

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

                let gene_id = ptir.source_geneid;

                let gene_idx = global_string_id.get_or_insert_gene(&gene_id);

                if (gene_idx as usize) == geneid_indexes.len() {
                    geneid_indexes.push(GeneIndex::new(ptir.junctions(), ptir.start, ptir.end))
                } else if (gene_idx as usize) < geneid_indexes.len() {
                    geneid_indexes[gene_idx as usize].add(ptir.junctions(), ptir.start, ptir.end);
                } else {
                    error!("Gene index is larger than current GeneIndex vector length.");
                    std::process::exit(1)
                }

                if let Some(juncs) = ptir.junctions() {
                    if matches!(ptir.strand, ISOMSTRAND::Plus) {
                        for j in juncs {
                            junction_plus.push(Junction {
                                start: j.0,
                                end: j.1,
                            })
                        }
                    } else {
                        for j in juncs {
                            junction_minus.push(Junction {
                                start: j.0,
                                end: j.1,
                            })
                        }
                    }
                };

                if ptir.n_exons == 1 {
                    interval_mono_exon.push(Interval {
                        start: ptir.start,
                        stop: ptir.end,
                        val: ptir_idx,
                    });
                } else {
                    interval_multi_exon.push(Interval {
                        start: ptir.start,
                        stop: ptir.end,
                        val: ptir_idx,
                    });
                }

                ptirs.push(RefPTIR {
                    base: ptir,
                    attrs: span_vec,
                    gene_name,
                });
            }

            chrom_indexes.push(ChromIndex::new(
                junction_plus,
                junction_minus,
                interval_mono_exon,
                interval_multi_exon,
            ))
        }

        attr_string_pool.shrink_to_read_only();

        Ok(Self {
            filename: gtf_path,
            ptirs,
            attr_string_pool,
            chroms: chrom_indexes,
            genes: geneid_indexes,
            stringids: global_string_id,
        })
    }

    pub fn find_ovlp_from_mono_refs(&self, chr_name: &str, start: u32, end: u32) -> Option<Vec<&RefPTIR>> {
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

    pub fn find_ovlp_from_multi_refs(&self, chr_name: &str, start: u32, end: u32) -> Option<Vec<&RefPTIR>> {
        todo!()
    }
}

impl LogMemSize for RefPTIRManager {
    fn get_mem_size(&self) -> usize {
        todo!()
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
