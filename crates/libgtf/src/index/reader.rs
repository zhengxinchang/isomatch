use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufReader, Cursor, ErrorKind, Read, Seek, SeekFrom},
    path::Path,
};

// use log::{error, warn};

use crate::traits::{Decodable, DiskSize, PartialLoad};

use crate::{
    error::{Error, Result},
    index::{
        junction_pool::JunctionPool, splice_site_pool::SpliceSitePool, string_pool::StringPool,
        tx::TxBase, tx::TxBaseLoadArgs,
    },
};
// use super::format::{ChromDirectoryEntry, IndexHeader};
use crate::index::{ChromDirectoryEntry, IndexHeader};

pub struct IndexReader {
    file_id: usize,
    header: IndexHeader,
    chroms: Vec<ChromDirectoryEntry>,
    /// Chrom names in chrom_id order (index = chrom_id - 1).
    chrom_names: Vec<String>,
    /// Map from chrom name to chrom_id for fast lookup.
    chrom_name_to_id: HashMap<String, u16>,
    /// Seqids present in the source GTF but absent from the reference FASTA.
    missing_seqids: Vec<String>,
    file: File,
}

impl IndexReader {
    pub fn load_header<P: AsRef<Path>>(path: P) -> Result<IndexHeader> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        IndexHeader::decode_from(&mut reader, ())
    }

    pub fn open(file: File, file_id: usize) -> Result<IndexReader> {
        let mut reader = BufReader::new(file);
        let header = IndexHeader::decode_from(&mut reader, ())?;

        // if header.version != ISOMX_VERSION {

        //     return Err(Error::UnsupportedVersion { found: header.version, expected: ISOMX_VERSION });
        // }

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
                    "invalid chrom_id 0 in directory entry",
                )
            })? as usize;

            if chrom_idx >= chrom_names.len() {
                return Err(Error::InvalidIndex {
                    reason: "Chromosome names does not match chromosome ids".to_string(),
                });
            }

            let start = entry.chrom_name_offset as usize;
            let end = start + entry.chrom_name_len as usize;
            if end > chrom_name_table.len() {
                return Err(Error::InvalidIndex {
                    reason: format!(
                        "chrom name slice [{}..{}) exceeds name table length {}",
                        start,
                        end,
                        chrom_name_table.len()
                    ),
                });
            }

            let chrom_name = std::str::from_utf8(&chrom_name_table[start..end])
                .map_err(|e| io::Error::new(ErrorKind::InvalidData, e.to_string()))?
                .to_string();

            if chrom_name_to_id
                .insert(chrom_name.clone(), entry.chrom_id)
                .is_some()
            {
                return Err(Error::InvalidIndex {
                    reason: format!("duplicate chromosome name in index: {}", chrom_name),
                });
            }

            chrom_names[chrom_idx] = chrom_name;
        }

        let mut missing_seqids = Vec::with_capacity(header.missing_seqid_count as usize);
        if header.missing_seqid_table_len > 0 {
            let mut table = vec![0u8; header.missing_seqid_table_len as usize];
            reader.read_exact(&mut table)?;
            let mut pos = 0usize;
            while pos + 2 <= table.len() {
                let len = u16::from_le_bytes(table[pos..pos + 2].try_into().unwrap()) as usize;
                pos += 2;
                if pos + len > table.len() {
                    return Err(Error::InvalidIndex {
                        reason: "missing seqid table entry exceeds table bounds".to_string(),
                    });
                }
                let name = std::str::from_utf8(&table[pos..pos + len])
                    .map_err(|e| io::Error::new(ErrorKind::InvalidData, e.to_string()))?
                    .to_string();
                missing_seqids.push(name);
                pos += len;
            }
            // warn!("Index skipped transcripts on missing reference seqid(s):",);
            // warn!("{}", missing_seqids.join(","));
            // warn!(
            //     "Those transcripts will not be processed. You may consider redo index step with proper reference genome FASTA."
            // );
        }

        Ok(Self {
            file_id,
            header,
            chroms,
            chrom_names,
            chrom_name_to_id,
            missing_seqids,
            file: reader.into_inner(),
        })
    }

    pub fn get_chromosome_reader(&mut self, chrom_name: &str) -> Result<ChromBlockReader> {
        let chrom_id =
            *self
                .chrom_name_to_id
                .get(chrom_name)
                .ok_or_else(|| Error::ChromosomeNotFound {
                    name: chrom_name.to_string(),
                })?;

        let entry = self
            .chroms
            .iter()
            .find(|entry| entry.chrom_id == chrom_id)
            .ok_or_else(|| Error::InvalidIndex {
                reason: format!("missing directory entry for chromosome id {}", chrom_id),
            })?;

        let chrom_name = self.chrom_names[(chrom_id - 1) as usize].clone();

        ChromBlockReader::new(
            self.file_id,
            chrom_id,
            chrom_name,
            entry.global_tx_count,
            entry.global_junction_pool_offset,
            entry.global_junction_pool_len as usize,
            entry.global_string_pool_offset,
            entry.global_string_pool_len,
            entry.global_splice_site_pool_offset,
            entry.global_splice_site_pool_len as usize,
            self.file.try_clone()?,
            entry.global_tx_offset,
            0,
        )
    }

    pub fn get_chromosome_readers_map(&mut self) -> Result<HashMap<String, ChromBlockReader>> {
        let mut readers = HashMap::default();
        for chr_name in self.chrom_names.clone() {
            let reader =
                self.get_chromosome_reader(&chr_name)
                    .map_err(|e| Error::InvalidIndex {
                        reason: format!(
                            "Can ont get chromosome level data from index. Reason {:?}",
                            e
                        ),
                    })?;
            readers.insert(chr_name.to_string(), reader);
        }
        Ok(readers)
    }

    /// Scan all chromosomes and build a transcript_id → tx_idx lookup map.
    /// O(n) one-time cost
    pub fn build_txid_index(&mut self) -> Result<HashMap<String, u64>> {
        let mut map = HashMap::new();
        for chrom_name in self.chrom_names.clone() {
            let mut cr =
                self.get_chromosome_reader(&chrom_name)
                    .map_err(|e| Error::InvalidIndex {
                        reason: e.to_string(),
                    })?;
            while let Some(tx) =
                ChromBlockReader::next_record(&mut cr).map_err(|e| Error::InvalidIndex {
                    reason: e.to_string(),
                })?
            {
                let tx_id = tx.source_tx_id(&cr.string_pool);
                map.insert(tx_id, tx.tx_idx());
            }
        }
        Ok(map)
    }

    pub fn version(&self) -> u32 {
        self.header.version
    }

    pub fn md5(&self) -> [u8; 16] {
        self.header.md5
    }

    pub fn transcript_count(&self) -> u64 {
        self.header.total_tx_n
    }

    pub fn chromosome_names(&self) -> &[String] {
        &self.chrom_names
    }

    pub fn contains_chromosome(&self, name: &str) -> bool {
        self.chrom_name_to_id.contains_key(name)
    }

    pub fn missing_seqids(&self) -> &[String] {
        &self.missing_seqids
    }
}

pub struct ChromBlockReader {
    file_id: usize,
    chrom_id: u16,
    chrom_name: String,
    tx_count: u64,
    junction_pool: JunctionPool,
    string_pool: StringPool,
    splice_site_pool: SpliceSitePool,
    file: File,
    tx_base_offset: u64,
    next_tx_idx: u64,
}

impl ChromBlockReader {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        file_id: usize,
        chrom_id: u16,
        chrom_name: String,
        tx_count: u64,
        junction_pool_offset: u64,
        junction_pool_len: usize,
        string_pool_offset: u64,
        string_pool_len: u64,
        splice_site_pool_offset: u64,
        splice_site_pool_len: usize,
        mut file: File,
        tx_base_offset: u64,
        next_tx_idx: u64,
    ) -> Result<ChromBlockReader> {
        let junction_pool = ChromBlockReader::load_junction_pool(
            &mut file,
            junction_pool_offset,
            junction_pool_len,
        )?;
        let string_pool =
            ChromBlockReader::load_string_pool(&mut file, string_pool_offset, string_pool_len)?;
        let splice_site_pool = ChromBlockReader::load_splice_site_pool(
            &mut file,
            splice_site_pool_offset,
            splice_site_pool_len,
        )?;
        Ok(Self {
            file_id,
            chrom_id,
            chrom_name,
            tx_count,
            junction_pool,
            string_pool,
            splice_site_pool,
            file,
            tx_base_offset,
            next_tx_idx,
        })
    }

    pub fn next_record(&mut self) -> Result<Option<TxBase>> {
        if self.next_tx_idx >= self.tx_count {
            return Ok(None);
        }

        let tx_offset = self.tx_base_offset + self.next_tx_idx * TxBase::DISK_SIZE as u64;
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

    pub fn file_id(&self) -> usize {
        self.file_id
    }

    pub fn chromosome_name(&self) -> &str {
        &self.chrom_name
    }

    pub fn junction_pool(&self) -> &JunctionPool {
        &self.junction_pool
    }

    pub fn string_pool(&self) -> &StringPool {
        &self.string_pool
    }

    pub fn splice_site_pool(&self) -> &SpliceSitePool {
        &self.splice_site_pool
    }

    fn decompress(file: &mut File, offset: u64, compressed_len: u64) -> Result<Vec<u8>> {
        let compressed_len = usize::try_from(compressed_len).map_err(|_| {
            io::Error::new(
                ErrorKind::InvalidData,
                format!("compressed pool length {compressed_len} exceeded usize"),
            )
        })?;
        let mut compressed = vec![0u8; compressed_len];
        file.seek(SeekFrom::Start(offset))?;
        file.read_exact(&mut compressed)?;
        zstd::decode_all(compressed.as_slice()).map_err(|e| Error::InvalidIndex {
            reason: format!("zstd decode failed for pool at offset {offset}: {e}"),
        })
    }

    pub(crate) fn load_junction_pool(
        file: &mut File,
        junction_pool_offset: u64,
        junction_pool_len: usize,
    ) -> Result<JunctionPool> {
        let decompressed = Self::decompress(file, junction_pool_offset, junction_pool_len as u64)?;
        let len = decompressed.len();
        JunctionPool::load_range(&mut Cursor::new(decompressed), 0, len, 0)
            .map_err(|e| Error::Io(io::Error::new(ErrorKind::InvalidData, e.to_string())))
    }

    pub(crate) fn load_string_pool(
        file: &mut File,
        string_pool_offset: u64,
        string_pool_len: u64,
    ) -> Result<StringPool> {
        let decompressed = Self::decompress(file, string_pool_offset, string_pool_len)?;
        let len = decompressed.len();
        StringPool::load_range(&mut Cursor::new(decompressed), 0, len, ())
            .map_err(|e| Error::Io(io::Error::new(ErrorKind::InvalidData, e.to_string())))
    }

    pub(crate) fn load_splice_site_pool(
        file: &mut File,
        splice_site_pool_offset: u64,
        splice_site_pool_len: usize,
    ) -> Result<SpliceSitePool> {
        let decompressed =
            Self::decompress(file, splice_site_pool_offset, splice_site_pool_len as u64)?;
        let len = decompressed.len();
        SpliceSitePool::load_range(&mut Cursor::new(decompressed), 0, len, ())
            .map_err(|e| Error::Io(io::Error::new(ErrorKind::InvalidData, e.to_string())))
    }
}
