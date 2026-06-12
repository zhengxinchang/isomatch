use std::fs::File;
use std::io::{self, BufWriter, Seek, SeekFrom, Write};

fn u16_from_usize(value: usize, label: &str) -> io::Result<u16> {
    u16::try_from(value).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{label} {value} exceeded u16"),
        )
    })
}

fn compress(raw: Vec<u8>) -> std::io::Result<Vec<u8>> {
    zstd::encode_all(raw.as_slice(), 3)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))
}

fn u64_from_usize(value: usize, label: &str) -> io::Result<u64> {
    u64::try_from(value).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{label} {value} exceeded u64"),
        )
    })
}

fn encode_compressed<T: crate::traits::Encodable<Error = crate::core::core_error::TxBaseError>>(
    pool: &T,
    writer: &mut impl Write,
) -> std::io::Result<usize> {
    let mut raw = Vec::new();
    pool.encode_to(&mut raw)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
    let compressed = compress(raw)?;
    let len = compressed.len();
    writer.write_all(&compressed)?;
    Ok(len)
}

use crate::core::tx_base::TxBase;
use crate::index::format::{ChromBlockBuilder, ChromDirectoryEntry, IndexHeader};
use crate::traits::{DiskSize, Encodable};

pub struct IndexBuilder {
    header: IndexHeader,
    entries: Vec<ChromDirectoryEntry>,
    /// Pre-computed (offset_in_table, len) for each chrom, indexed by chrom_id - 1.
    chrom_name_offsets: Vec<(u32, u32)>,
    current_offset: u64,
    // total_tx_n: u32,
    file: BufWriter<File>,
}

impl IndexBuilder {
    pub fn new(
        file: File,
        chrom_names: Vec<String>,
        gtf_file_size: u64,
        md5: [u8; 16],
        has_ref_hash: bool,
        has_seq_hash: bool,
        missing_seqids: Vec<String>,
    ) -> std::io::Result<Self> {
        let mut file = BufWriter::new(file);
        let chrom_count = u32::try_from(chrom_names.len()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("chromosome count {} exceeded u32", chrom_names.len()),
            )
        })?;

        // Build the chrom name table bytes and pre-compute per-chrom offsets.
        let mut name_table: Vec<u8> = Vec::new();
        let mut chrom_name_offsets: Vec<(u32, u32)> = Vec::with_capacity(chrom_names.len());
        for name in &chrom_names {
            let offset = u32::try_from(name_table.len()).map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("chrom name table offset {} exceeded u32", name_table.len()),
                )
            })?;
            let len = u32::try_from(name.len()).map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("chrom name length {} exceeded u32", name.len()),
                )
            })?;
            name_table.extend_from_slice(name.as_bytes());
            chrom_name_offsets.push((offset, len));
        }
        let chrom_name_table_len = u32::try_from(name_table.len()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("chrom name table length {} exceeded u32", name_table.len()),
            )
        })?;

        // Build the missing seqid table: each entry is u16 len + utf-8 bytes.
        let mut missing_seqid_table: Vec<u8> = Vec::new();
        for name in &missing_seqids {
            let bytes = name.as_bytes();
            let len = u16_from_usize(bytes.len(), "missing seqid length")?;
            missing_seqid_table.extend_from_slice(&len.to_le_bytes());
            missing_seqid_table.extend_from_slice(bytes);
        }
        let missing_seqid_count = u32::try_from(missing_seqids.len()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("missing seqid count {} exceeded u32", missing_seqids.len()),
            )
        })?;
        let missing_seqid_table_len = u32::try_from(missing_seqid_table.len()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "missing seqid table length {} exceeded u32",
                    missing_seqid_table.len()
                ),
            )
        })?;

        let header = IndexHeader::new(
            chrom_count,
            gtf_file_size,
            0,
            md5,
            has_ref_hash,
            has_seq_hash,
            chrom_name_table_len,
            missing_seqid_count,
            missing_seqid_table_len,
        );

        // Write placeholder header (4 KB)
        file.write_all(&[0u8; IndexHeader::DISK_SIZE])?;
        // Write placeholder directory (N × DISK_SIZE B)
        file.write_all(&vec![
            0u8;
            chrom_count as usize * ChromDirectoryEntry::DISK_SIZE
        ])?;
        // Write chrom name table — fixed, never rewritten
        file.write_all(&name_table)?;
        // Write missing seqid table — fixed, never rewritten
        file.write_all(&missing_seqid_table)?;

        let current_offset = (IndexHeader::DISK_SIZE
            + chrom_count as usize * ChromDirectoryEntry::DISK_SIZE
            + name_table.len()
            + missing_seqid_table.len()) as u64;

        Ok(Self {
            header,
            entries: Vec::with_capacity(chrom_count as usize),
            chrom_name_offsets,
            current_offset,
            // total_tx_n: 0,
            file,
        })
    }

    pub fn add_chrom(&mut self, mut entry: ChromBlockBuilder) -> std::io::Result<()> {
        entry.finalize();

        let (chrom_name_offset, chrom_name_len) =
            self.chrom_name_offsets[(entry.chrom_id - 1) as usize];

        let tx_offset = self.current_offset;
        let tx_bytes = entry.txs.len() * TxBase::DISK_SIZE;

        for tx in &entry.txs {
            tx.encode_to(&mut self.file)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        }
        self.current_offset += tx_bytes as u64;

        let junction_pool_offset = tx_offset + tx_bytes as u64;
        let junction_pool_len = encode_compressed(&entry.junction_pool, &mut self.file)?;
        self.current_offset += u64_from_usize(junction_pool_len, "junction pool length")?;

    
        let string_pool_offset =
            junction_pool_offset + u64_from_usize(junction_pool_len, "junction pool length")?;
        let string_pool_len = encode_compressed(&entry.string_pool, &mut self.file)?;
        self.current_offset += u64_from_usize(string_pool_len, "string pool length")?;


        let splice_site_pool_offset =
            string_pool_offset + u64_from_usize(string_pool_len, "string pool length")?;
        let splice_site_pool_len = encode_compressed(&entry.splice_site_pool, &mut self.file)?;
        self.current_offset += u64_from_usize(splice_site_pool_len, "splice site pool length")?;


        let tx_count = u64::try_from(entry.txs.len()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "chromosome transcript count {} exceeded u64",
                    entry.txs.len()
                ),
            )
        })?;
        self.header.total_tx_n = self
            .header
            .total_tx_n
            .checked_add(tx_count)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "total transcript count exceeded u64",
                )
            })?;
        self.entries.push(ChromDirectoryEntry {
            chrom_id: entry.chrom_id,
            chrom_name_offset,
            chrom_name_len,
            global_tx_count: tx_count,
            global_tx_offset: tx_offset,
            global_junction_pool_offset: junction_pool_offset,
            global_junction_count: u64_from_usize(junction_pool_len, "junction pool length")?,
            global_string_pool_offset: string_pool_offset,
            global_string_len: u64_from_usize(string_pool_len, "string pool length")?,
            global_splice_site_pool_offset: splice_site_pool_offset,
            global_splice_site_pool_len: u64_from_usize(
                splice_site_pool_len,
                "splice site pool length",
            )?,
        });

        Ok(())
    }

    /// Seek back and write the real header and directory.
    pub fn finalize(mut self) -> std::io::Result<()> {
        self.header.index_file_size = self.current_offset;
        // self.header.total_tx_n = self.total_tx_n;
        self.file.seek(SeekFrom::Start(0))?;
        self.header.encode_to(&mut self.file)?;

        self.file
            .seek(SeekFrom::Start(IndexHeader::DISK_SIZE as u64))?;
        for entry in &self.entries {
            entry.encode_to(&mut self.file)?;
        }

        self.file.flush()
    }
}
