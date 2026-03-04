use std::{
    collections::HashMap,
    fs::File,
    io::{BufReader, Read},
};

use crate::index::format::{ChromDirectoryEntry, IndexHeader};
use crate::traits::Decodable;

pub struct Index {
    pub header: IndexHeader,
    pub chroms: Vec<ChromDirectoryEntry>,
    /// Chrom names in chrom_id order (index = chrom_id - 1).
    pub chrom_names: Vec<String>,
    /// Map from chrom name to chrom_id for fast lookup.
    pub chrom_name_to_id: HashMap<String, u16>,
    pub file: File,
}

impl Index {
    pub fn open(file: File) -> Result<Self, std::io::Error> {
        let mut reader = BufReader::new(file);

        let header = IndexHeader::decode_from(&mut reader, ())?;

        let mut chroms = Vec::with_capacity(header.chrom_count as usize);
        for _ in 0..header.chrom_count {
            chroms.push(ChromDirectoryEntry::decode_from(&mut reader, ())?);
        }

        // Chrom name table is contiguous right after the directory — one sequential read.
        let mut name_table = vec![0u8; header.chrom_name_table_len as usize];
        reader.read_exact(&mut name_table)?;

        let mut chrom_names = Vec::with_capacity(header.chrom_count as usize);
        let mut chrom_name_to_id = HashMap::with_capacity(header.chrom_count as usize);
        for entry in &chroms {
            let start = entry.chrom_name_offset as usize;
            let end = start + entry.chrom_name_len as usize;
            let name = std::str::from_utf8(&name_table[start..end])
                .map_err(|_| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "invalid UTF-8 in chrom name",
                    )
                })?
                .to_string();
            chrom_name_to_id.insert(name.clone(), entry.chrom_id);
            chrom_names.push(name);
        }

        Ok(Self {
            header,
            chroms,
            chrom_names,
            chrom_name_to_id,
            file: reader.into_inner(),
        })
    }
}
