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
    core::{ptir::PTIR, string_pool::StringPool, tx_strand::ISOMSTRAND},
    index::reader::IndexReader,
    traits::LogMemSize,
};

type ChromId = u32;
type GeneId = u32;
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
}
struct ChromIndex {
    splice_site_starts: Vec<Pos>,
    splice_site_ends: Vec<Pos>,

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

        let mut max_end = 0;
        for junction in junctions_plus.iter() {
            starts.insert(junction.start);
            ends.insert(junction.end);
            max_end = max_end.max(junction.end);
            prefix_max_starts_plus.push(max_end);
        }

        let mut max_end = 0;
        for junction in junctions_minus.iter() {
            starts.insert(junction.start);
            ends.insert(junction.end);
            max_end = max_end.max(junction.end);
            prefix_max_starts_minus.push(max_end);
        }

        let mut splice_site_starts: Vec<Pos> = starts.into_iter().collect();
        splice_site_starts.sort();

        let mut splice_site_ends: Vec<Pos> = ends.into_iter().collect();
        splice_site_ends.sort();

        Self {
            splice_site_starts,
            splice_site_ends,
            junctions_plus,
            junctions_minus,
            prefix_max_starts_plus,
            prefix_max_starts_minus,
            refs_monoexon: Lapper::new(interval_mono_exon),
            refs_multiexon: Lapper::new(interval_multi_exon),
        }
    }

    pub fn has_junction(&self, junction: &Junction, strand: &ISOMSTRAND) -> bool {
        match strand {
            ISOMSTRAND::Minus => match self.junctions_minus.binary_search(junction) {
                Ok(_) => return true,
                Err(_) => return false,
            },
            ISOMSTRAND::Plus => match self.junctions_plus.binary_search(junction) {
                Ok(_) => return true,
                Err(_) => return false,
            },
            ISOMSTRAND::Unknown => {
                panic!("This should not happen as it only acccept stranded transcript in new().");
            }
        }
    }

    pub fn has_splice_site(&self, start: u32, end: u32) -> (bool, bool) {
        let has_start = self.splice_site_starts.binary_search(&start);
        let has_end = self.splice_site_ends.binary_search(&end);
        (has_start.is_ok(), has_end.is_ok())
    }

    pub fn with_in_junction(&self, txboundary: &Junction, strand: &ISOMSTRAND) -> bool {
        let (junctions, prefix_max_ends) = match strand {
            ISOMSTRAND::Plus => (&self.junctions_plus, &self.prefix_max_starts_plus),
            ISOMSTRAND::Minus => (&self.junctions_minus, &self.prefix_max_starts_minus),
            ISOMSTRAND::Unknown => {
                panic!("This should not happen as it only acccept stranded transcript in new().");
            }
        };

        let candidate_end = junctions.partition_point(|junction| junction.start <= txboundary.end);
        if candidate_end == 0 {
            return false;
        }

        let candidate_start =
            prefix_max_ends[..candidate_end].partition_point(|&max_end| max_end < txboundary.start);

        junctions[candidate_start..candidate_end]
            .iter()
            .any(|junction| junction.start <= txboundary.start && txboundary.end <= junction.end)
    }

    pub fn junction_inside_boundary(&self, txboundary: &Junction, strand: &ISOMSTRAND) -> bool {
        let junctions = match strand {
            ISOMSTRAND::Plus => &self.junctions_plus,
            ISOMSTRAND::Minus => &self.junctions_minus,
            ISOMSTRAND::Unknown => {
                panic!("This should not happen as it only acccept stranded transcript in new().");
            }
        };

        let candidate_start =
            junctions.partition_point(|junction| junction.start < txboundary.start);
        let candidate_end = junctions.partition_point(|junction| junction.start <= txboundary.end);

        junctions[candidate_start..candidate_end]
            .iter()
            .any(|junction| {
                txboundary.start <= junction.start
                    && junction.start < junction.end
                    && junction.end < txboundary.end
            })
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
        self.ends.push(end);
    }

    fn finalize(&mut self) {
        self.junctions.sort();
        self.junctions.dedup();
        self.starts.sort();
        self.starts.dedup();
        self.ends.sort();
        self.ends.dedup();
    }

    fn junction_pairs(&self) -> Vec<(u32, u32)> {
        self.junctions
            .iter()
            .map(|junction| (junction.start, junction.end))
            .collect()
    }

    fn starts(&self) -> &[Pos] {
        &self.starts
    }

    fn ends(&self) -> &[Pos] {
        &self.ends
    }
}
pub struct RefPTIRManager {
    pub filename: PathBuf,
    pub ptirs: Vec<RefPTIR>,
    pub attr_string_pool: StringPool,
    chroms: Vec<ChromIndex>,
    genes: Vec<GeneIndex>,
    stringids: StringID,
}

impl RefPTIRManager {
    pub fn open<P: AsRef<Path>>(gtf_path: P) -> Result<Self, ClassifyError> {
        let gtf_path = gtf_path.as_ref().to_path_buf();
        let mut isomx_path = gtf_path.clone();
        isomx_path.add_extension("isomx");

        let f = File::open(&isomx_path).map_err(|e| ClassifyError::FailedParseGTF {
            reason: format!(
                "Can not read reference GTF {}, reason: {}",
                &isomx_path.display(),
                e.to_string()
            ),
        })?;

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
            global_string_id.get_or_insert_chrom(chr_name);
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
                        "Reference GTF contains unstranded transcript: {}, not used in classify.",
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

                let gene_id = &ptir.source_geneid;

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
                    chr_name: chr_name.clone(),
                });
            }

            chrom_indexes.push(ChromIndex::new(
                junction_plus,
                junction_minus,
                interval_mono_exon,
                interval_multi_exon,
            ))
        }

        // Gene-level vectors are appended transcript-by-transcript while reading
        // the index. Finalize once so downstream classification sees deduplicated
        // known junction/start/end catalogs.
        for gene in &mut geneid_indexes {
            gene.finalize();
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

    pub fn has_chr(&self, chr_name: &str) -> bool {
        match self.stringids.chrom_id(chr_name) {
            Some(_) => true,
            None => false,
        }
    }

    pub fn find_ovlp_from_mono_refs(
        &self,
        chr_name: &str,
        start: u32,
        end: u32,
    ) -> Option<Vec<&RefPTIR>> {
        let chrom_id = self.stringids.chrom_id(chr_name)? as usize;
        let chrom = self.chroms.get(chrom_id)?;
        let results: Vec<&RefPTIR> = chrom
            .refs_monoexon
            .find(start, end)
            .map(|iv| &self.ptirs[iv.val])
            .collect();
        if results.is_empty() {
            None
        } else {
            Some(results)
        }
    }

    pub fn find_ovlp_from_multi_refs(
        &self,
        chr_name: &str,
        start: u32,
        end: u32,
    ) -> Option<Vec<&RefPTIR>> {
        let chrom_id = self.stringids.chrom_id(chr_name)? as usize;
        let chrom = self.chroms.get(chrom_id)?;
        let results: Vec<&RefPTIR> = chrom
            .refs_multiexon
            .find(start, end)
            .map(|iv| &self.ptirs[iv.val])
            .collect();
        if results.is_empty() {
            None
        } else {
            Some(results)
        }
    }

    pub fn find_overlapping_refs(&self, chr_name: &str, start: u32, end: u32) -> Vec<&RefPTIR> {
        let Some(chrom_id) = self.stringids.chrom_id(chr_name) else {
            return Vec::new();
        };
        let Some(chrom) = self.chroms.get(chrom_id as usize) else {
            return Vec::new();
        };

        let mut results: Vec<&RefPTIR> = chrom
            .refs_multiexon
            .find(start, end)
            .map(|iv| &self.ptirs[iv.val])
            .collect();
        results.extend(
            chrom
                .refs_monoexon
                .find(start, end)
                .map(|iv| &self.ptirs[iv.val]),
        );
        results
    }

    pub fn gene_junctions(&self, gene_id: &str) -> Option<Vec<(u32, u32)>> {
        let gene_idx = self.stringids.gene_id(gene_id)? as usize;
        self.genes.get(gene_idx).map(GeneIndex::junction_pairs)
    }

    pub fn gene_starts_ends(&self, gene_id: &str) -> Option<(&[u32], &[u32])> {
        let gene_idx = self.stringids.gene_id(gene_id)? as usize;
        let gene = self.genes.get(gene_idx)?;
        Some((gene.starts(), gene.ends()))
    }

    pub fn junction_match(
        &self,
        chr_name: &str,
        junctions: &[(u32, u32)],
        strand: &ISOMSTRAND,
    ) -> (bool, Vec<bool>) {
        let chrid = self.stringids.chrom_id(chr_name).unwrap();

        let mut is_all_junction_found = true;
        let mut junction_found_vec = Vec::new();
        for junction in junctions {
            let j = Junction {
                start: junction.0,
                end: junction.1,
            };
            let found = self.chroms[chrid as usize].has_junction(&j, strand);
            is_all_junction_found &= found;
            junction_found_vec.push(found);
        }

        (is_all_junction_found, junction_found_vec)
    }

    pub fn splice_site_match(
        &self,
        chr_name: &str,
        junctions: &[(u32, u32)],
    ) -> (bool, Vec<bool>, Vec<bool>) {
        let chrid = self.stringids.chrom_id(chr_name).unwrap();

        let mut is_all_splice_sites_found = true;
        let mut left_site_found_vec = Vec::with_capacity(junctions.len());
        let mut right_site_found_vec = Vec::with_capacity(junctions.len());

        for junction in junctions {
            let (left_found, right_found) =
                self.chroms[chrid as usize].has_splice_site(junction.0, junction.1);
            is_all_splice_sites_found &= left_found && right_found;
            left_site_found_vec.push(left_found);
            right_site_found_vec.push(right_found);
        }

        (
            is_all_splice_sites_found,
            left_site_found_vec,
            right_site_found_vec,
        )
    }

    pub fn contained_in_known_intron(
        &self,
        chr_name: &str,
        strand: &ISOMSTRAND,
        start: u32,
        end: u32,
    ) -> bool {
        let junction = Junction { start, end };
        let chrid = self.stringids.chrom_id(chr_name).unwrap();
        self.chroms[chrid as usize].with_in_junction(&junction, strand)
    }

    pub fn has_intron_retention_against_catalog(
        &self,
        chr_name: &str,
        strand: &ISOMSTRAND,
        exons: &[(u32, u32)],
    ) -> bool {
        let chrid = self.stringids.chrom_id(chr_name).unwrap();
        let chrom = &self.chroms[chrid as usize];

        exons.iter().any(|&(start, end)| {
            let exon = Junction { start, end };
            chrom.junction_inside_boundary(&exon, strand)
        })
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
