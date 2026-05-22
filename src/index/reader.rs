use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufReader, Cursor, ErrorKind, Read, Seek, SeekFrom},
};

use crate::traits::{Decodable, DiskSize, PartialLoad};
use crate::{
    core::{
        junction_pool::JunctionPool, splice_site_pool::SpliceSitePool, string_pool::StringPool,
        tx_base::TxBase, tx_base_impl::TxBaseLoadArgs,
    },
    index::index_error::IndexError,
};

use super::format::{ChromDirectoryEntry, IndexHeader};

pub struct IndexReader {
    pub file_id: usize,
    pub header: IndexHeader,
    pub chroms: Vec<ChromDirectoryEntry>,
    /// Chrom names in chrom_id order (index = chrom_id - 1).
    pub chrom_names: Vec<String>,
    /// Map from chrom name to chrom_id for fast lookup.
    pub chrom_name_to_id: HashMap<String, u16>,
    pub file: File,
}

impl IndexReader {
    pub fn open(file: File, file_id: usize) -> io::Result<Self> {
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
            file_id,
            header,
            chroms,
            chrom_names,
            chrom_name_to_id,
            file: reader.into_inner(),
        })
    }

    pub fn get_chromosome_reader(&mut self, chrom_name: &str) -> io::Result<ChromBlockReader> {
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

        ChromBlockReader::new(
            self.file_id,
            chrom_id,
            chrom_name,
            entry.global_tx_count,
            entry.global_junction_pool_offset as u64,
            entry.global_junction_count as usize,
            entry.global_string_pool_offset as u64,
            entry.global_string_len as usize,
            entry.global_splice_site_pool_offset as u64,
            entry.global_splice_site_pool_len as usize,
            self.file.try_clone()?,
            entry.global_tx_offset as u64,
            0,
        )
    }

    pub fn get_chromosome_readers_map(
        &mut self,
    ) -> Result<HashMap<String, ChromBlockReader>, IndexError> {
        let mut readers = HashMap::default();
        for chr_name in self.chrom_names.clone() {
            let reader =
                self.get_chromosome_reader(&chr_name)
                    .map_err(|e| IndexError::FailReadIndex {
                        reason: format!(
                            "Can ont get chromosome level data from index. Reason {:?}",
                            e
                        ),
                    })?;
            readers.insert(chr_name.to_string(), reader);
        }
        Ok(readers)
    }
}

pub struct ChromBlockReader {
    pub file_id: usize,
    pub chrom_id: u16,
    pub chrom_name: String,
    pub tx_count: u32,
    pub junction_pool_offset: u64,
    pub junction_pool_len: usize,
    pub junction_pool: JunctionPool,
    pub string_pool_offset: u64,
    pub string_pool_len: usize,
    pub string_pool: StringPool,
    pub splice_site_pool_offset: u64,
    pub splice_site_pool_len: usize,
    pub splice_site_pool: SpliceSitePool,
    file: File,
    tx_base_offset: u64,
    next_tx_idx: u32,
}

impl ChromBlockReader {
    pub fn new(
        file_id: usize,
        chrom_id: u16,
        chrom_name: String,
        tx_count: u32,
        junction_pool_offset: u64,
        junction_pool_len: usize,
        string_pool_offset: u64,
        string_pool_len: usize,
        splice_site_pool_offset: u64,
        splice_site_pool_len: usize,
        mut file: File,
        tx_base_offset: u64,
        next_tx_idx: u32,
    ) -> io::Result<ChromBlockReader> {
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
            junction_pool_offset,
            junction_pool_len,
            junction_pool: junction_pool,
            string_pool_offset,
            string_pool_len,
            string_pool: string_pool,
            splice_site_pool_offset,
            splice_site_pool_len,
            splice_site_pool: splice_site_pool,
            file,
            tx_base_offset,
            next_tx_idx,
        })
    }

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

    fn decompress(file: &mut File, offset: u64, compressed_len: usize) -> io::Result<Vec<u8>> {
        let mut compressed = vec![0u8; compressed_len];
        file.seek(SeekFrom::Start(offset))?;
        file.read_exact(&mut compressed)?;
        zstd::decode_all(compressed.as_slice())
            .map_err(|e| io::Error::new(ErrorKind::InvalidData, e.to_string()))
    }

    pub fn load_junction_pool(
        file: &mut File,
        junction_pool_offset: u64,
        junction_pool_len: usize,
    ) -> io::Result<JunctionPool> {
        let decompressed = Self::decompress(file, junction_pool_offset, junction_pool_len)?;
        let len = decompressed.len();
        JunctionPool::load_range(&mut Cursor::new(decompressed), 0, len, 0)
            .map_err(|e| io::Error::new(ErrorKind::InvalidData, e.to_string()))
    }

    pub fn load_string_pool(
        file: &mut File,
        string_pool_offset: u64,
        string_pool_len: usize,
    ) -> io::Result<StringPool> {
        let decompressed = Self::decompress(file, string_pool_offset, string_pool_len)?;
        let len = decompressed.len();
        StringPool::load_range(&mut Cursor::new(decompressed), 0, len, ())
            .map_err(|e| io::Error::new(ErrorKind::InvalidData, e.to_string()))
    }

    pub fn load_splice_site_pool(
        file: &mut File,
        splice_site_pool_offset: u64,
        splice_site_pool_len: usize,
    ) -> io::Result<SpliceSitePool> {
        let decompressed = Self::decompress(file, splice_site_pool_offset, splice_site_pool_len)?;
        let len = decompressed.len();
        SpliceSitePool::load_range(&mut Cursor::new(decompressed), 0, len, ())
            .map_err(|e| io::Error::new(ErrorKind::InvalidData, e.to_string()))
    }
}

impl Iterator for ChromBlockReader {
    type Item = TxBase;

    fn next(&mut self) -> Option<Self::Item> {
        match ChromBlockReader::next(self) {
            Ok(txbase) => txbase,
            Err(e) => {
                eprintln!("cannot read next transcript from isomx index: {}", e);
                std::process::exit(1);
            }
        }
    }
}
