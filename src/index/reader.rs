use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufReader, ErrorKind, Read},
};

use crate::core::{
    junction_pool::JunctionPool, splice_site_pool::SpliceSitePool, string_pool::StringPool,
    tx_base::TxBase, tx_base_impl::TxBaseLoadArgs,
};
use crate::traits::{Decodable, DiskSize, PartialLoad};

use super::format::{ChromDirectoryEntry, IndexHeader};

pub struct IndexReader {
    pub header: IndexHeader,
    pub chroms: Vec<ChromDirectoryEntry>,
    /// Chrom names in chrom_id order (index = chrom_id - 1).
    pub chrom_names: Vec<String>,
    /// Map from chrom name to chrom_id for fast lookup.
    pub chrom_name_to_id: HashMap<String, u16>,
    pub file: File,
}

pub struct ChromBlockReader {
    pub chrom_id: u16,
    pub chrom_name: String,
    pub tx_count: u32,
    pub junction_pool_offset: u64,
    pub junction_pool_len: usize,
    pub string_pool_offset: u64,
    pub string_pool_len: usize,
    pub splice_site_pool_offset: u64,
    pub splice_site_pool_len: usize,
    file: File,
    tx_base_offset: u64,
    next_tx_idx: u32,
}

impl IndexReader {
    pub fn open(file: File) -> io::Result<Self> {
        let mut reader = BufReader::new(file);
        let header = IndexHeader::decode_from(&mut reader, ())?;

        let mut chroms = Vec::with_capacity(header.chrom_count as usize);
        for _ in 0..header.chrom_count {
            chroms.push(ChromDirectoryEntry::decode_from(&mut reader, ())?);
        }

        let mut chrom_name_table = vec![0u8; header.chrom_name_table_len as usize];
        reader.read_exact(&mut chrom_name_table)?;

        let mut chrom_names = vec![String::new(); header.chrom_count as usize];
        let mut chrom_name_to_id = HashMap::with_capacity(header.chrom_count as usize);

        for entry in &chroms {
            let chrom_idx = entry.chrom_id.checked_sub(1).ok_or_else(|| {
                io::Error::new(
                    ErrorKind::InvalidData,
                    format!("invalid chrom_id 0 in directory entry"),
                )
            })? as usize;

            if chrom_idx >= chrom_names.len() {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    format!(
                        "chrom_id {} exceeds declared chrom_count {}",
                        entry.chrom_id, header.chrom_count
                    ),
                ));
            }

            let start = entry.chrom_name_offset as usize;
            let end = start + entry.chrom_name_len as usize;
            if end > chrom_name_table.len() {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    format!(
                        "chrom name slice [{}..{}) exceeds name table length {}",
                        start,
                        end,
                        chrom_name_table.len()
                    ),
                ));
            }

            let chrom_name = std::str::from_utf8(&chrom_name_table[start..end])
                .map_err(|e| io::Error::new(ErrorKind::InvalidData, e.to_string()))?
                .to_string();

            if chrom_name_to_id
                .insert(chrom_name.clone(), entry.chrom_id)
                .is_some()
            {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    format!("duplicate chromosome name in index: {}", chrom_name),
                ));
            }

            chrom_names[chrom_idx] = chrom_name;
        }

        Ok(Self {
            header,
            chroms,
            chrom_names,
            chrom_name_to_id,
            file: reader.into_inner(),
        })
    }

    pub fn get_chrosome_reader(&mut self, chrom_name: &str) -> io::Result<ChromBlockReader> {
        let chrom_id = *self.chrom_name_to_id.get(chrom_name).ok_or_else(|| {
            io::Error::new(
                ErrorKind::NotFound,
                format!("chromosome not found in index: {}", chrom_name),
            )
        })?;

        let entry = self
            .chroms
            .iter()
            .find(|entry| entry.chrom_id == chrom_id)
            .ok_or_else(|| {
                io::Error::new(
                    ErrorKind::InvalidData,
                    format!("missing directory entry for chromosome id {}", chrom_id),
                )
            })?;

        let chrom_name = self.chrom_names[(chrom_id - 1) as usize].clone();

        Ok(ChromBlockReader {
            chrom_id,
            chrom_name,
            tx_count: entry.global_tx_count,
            junction_pool_offset: entry.global_junction_pool_offset as u64,
            junction_pool_len: entry.global_junction_count as usize,
            string_pool_offset: entry.global_string_pool_offset as u64,
            string_pool_len: entry.global_string_len as usize,
            splice_site_pool_offset: entry.global_splice_site_pool_offset as u64,
            splice_site_pool_len: entry.global_splice_site_pool_len as usize,
            file: self.file.try_clone()?,
            tx_base_offset: entry.global_tx_offset as u64,
            next_tx_idx: 0,
        })
    }

    pub fn get_chromosome_reader(&mut self, chrom_name: &str) -> io::Result<ChromBlockReader> {
        self.get_chrosome_reader(chrom_name)
    }
}

impl ChromBlockReader {
    pub fn next(&mut self) -> io::Result<Option<TxBase>> {
        if self.next_tx_idx >= self.tx_count {
            return Ok(None);
        }

        let tx_offset = self.tx_base_offset + self.next_tx_idx as u64 * TxBase::DISK_SIZE as u64;
        let tx = TxBase::load_range(
            &mut self.file,
            tx_offset,
            TxBase::DISK_SIZE,
            TxBaseLoadArgs {
                chrom_id: self.chrom_id,
            },
        )
        .map_err(|e| io::Error::new(ErrorKind::InvalidData, e.to_string()))?;

        self.next_tx_idx += 1;
        Ok(Some(tx))
    }

    pub fn reset(&mut self) {
        self.next_tx_idx = 0;
    }

    pub fn load_junction_pool(&mut self) -> io::Result<JunctionPool> {
        JunctionPool::load_range(
            &mut self.file,
            self.junction_pool_offset,
            self.junction_pool_len,
            0,
        )
        .map_err(|e| io::Error::new(ErrorKind::InvalidData, e.to_string()))
    }

    pub fn load_string_pool(&mut self) -> io::Result<StringPool> {
        StringPool::load_range(
            &mut self.file,
            self.string_pool_offset,
            self.string_pool_len,
            (),
        )
        .map_err(|e| io::Error::new(ErrorKind::InvalidData, e.to_string()))
    }

    pub fn load_splice_site_pool(&mut self) -> io::Result<SpliceSitePool> {
        SpliceSitePool::load_range(
            &mut self.file,
            self.splice_site_pool_offset,
            self.splice_site_pool_len,
            (),
        )
        .map_err(|e| io::Error::new(ErrorKind::InvalidData, e.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn verify_index_readback() {
        let path = Path::new("test/isoseq_transcripts.sorted.filtered_lite.clean.isomx");
        if !path.exists() {
            eprintln!("Skipping test: index file not found");
            return;
        }

        let file = File::open(path).expect("cannot open index file");
        let mut reader = IndexReader::open(file).expect("cannot open index");

        assert_eq!(&reader.header.magic, b"ISOM");
        assert_eq!(reader.header.version, 1);
        assert!(reader.header.chrom_count > 0);
        assert_eq!(reader.chrom_names.len(), reader.header.chrom_count as usize);

        let mut total_tx = 0u32;
        let mut total_canonical = 0u32;
        let mut total_splice_sites = 0u32;

        for chrom_name in reader.chrom_names.clone() {
            let mut cr = reader
                .get_chromosome_reader(&chrom_name)
                .expect("cannot get chrom reader");

            let jp = cr.load_junction_pool().expect("cannot load junction pool");
            let sp = cr.load_string_pool().expect("cannot load string pool");
            let ssp = cr
                .load_splice_site_pool()
                .expect("cannot load splice site pool");

            while let Some(tx) = cr.next().expect("cannot read tx") {
                total_tx += 1;

                // validate string spans
                let tx_id = sp.get(tx.tx_id_span).expect("bad tx_id span");
                let gene_id = sp.get(tx.gene_id_span).expect("bad gene_id span");
                assert!(!tx_id.is_empty());
                assert!(!gene_id.is_empty());

                // validate junction span
                let junctions = jp.get(tx.junctions_span).expect("bad junction span");
                assert_eq!(junctions.len(), tx.junctions_span.count as usize);

                // validate splice site span
                if !tx.splice_sites_span.is_empty() {
                    let sites = ssp
                        .get_pair(tx.splice_sites_span)
                        .expect("bad splice site span");
                    assert_eq!(sites.len(), tx.splice_sites_span.count as usize);
                    // splice sites count should match junction pairs (n_exons - 1)
                    assert_eq!(sites.len(), (tx.n_exons - 1) as usize);

                    for site in sites {
                        total_splice_sites += 1;
                        if site.is_canonical() {
                            total_canonical += 1;
                        }
                    }
                } else {
                    // mono-exon transcripts have no splice sites
                    assert_eq!(tx.n_exons, 1);
                }
            }
        }

        assert!(total_tx > 0, "no transcripts found");
        eprintln!("Total transcripts: {}", total_tx);
        eprintln!("Total splice sites: {}", total_splice_sites);
        eprintln!("Canonical (GT-AG): {}", total_canonical);
        if total_splice_sites > 0 {
            let rate = total_canonical as f64 / total_splice_sites as f64 * 100.0;
            eprintln!("Canonical rate: {:.2}%", rate);
            // expect a high canonical rate for real data
            assert!(rate > 90.0, "canonical rate suspiciously low: {:.2}%", rate);
        }
    }
}
