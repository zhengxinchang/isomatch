use std::fs::File;
use std::io::{self, BufWriter, Seek, SeekFrom, Write};

fn u32_from_usize(value: usize, label: &str) -> io::Result<u32> {
    u32::try_from(value).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{label} {value} exceeded u32"),
        )
    })
}

fn u32_from_u64(value: u64, label: &str) -> io::Result<u32> {
    u32::try_from(value).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{label} {value} exceeded u32"),
        )
    })
}

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

fn encode_compressed<T: crate::traits::Encodable<Error = crate::core::core_error::TxBaseError>>(
    pool: &T,
    writer: &mut impl Write,
) -> std::io::Result<u32> {
    let mut raw = Vec::new();
    pool.encode_to(&mut raw)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
    let compressed = compress(raw)?;
    let len = u32_from_usize(compressed.len(), "compressed pool length")?;
    writer.write_all(&compressed)?;
    Ok(len)
}

use crate::core::tx_base::TxBase;
use crate::index::format::{ChromBlockBuilder, ChromDirectoryEntry, IndexHeader};
use crate::traits::{DiskSize, Encodable};

/// 1. new() write 4k header place holder,  + N×34B  Directory + Chrom Name Table + Missing SeqID Table
/// 2. add_chrom() write per-chrom data. drop immediately after it is done.
/// 3. finalize() seek back header place holder, fill up with Header and Directory
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
    /// `chrom_names` must be in the same order as chroms will appear in the GTF
    /// (i.e. the order returned by `profile_gtf`).
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
        let chrom_count = u32_from_usize(chrom_names.len(), "chromosome count")?;

        // Build the chrom name table bytes and pre-compute per-chrom offsets.
        let mut name_table: Vec<u8> = Vec::new();
        let mut chrom_name_offsets: Vec<(u32, u32)> = Vec::with_capacity(chrom_names.len());
        for name in &chrom_names {
            let offset = u32_from_usize(name_table.len(), "chrom name table offset")?;
            let len = u32_from_usize(name.len(), "chrom name length")?;
            name_table.extend_from_slice(name.as_bytes());
            chrom_name_offsets.push((offset, len));
        }
        let chrom_name_table_len = u32_from_usize(name_table.len(), "chrom name table length")?;

        // Build the missing seqid table: each entry is u16 len + utf-8 bytes.
        let mut missing_seqid_table: Vec<u8> = Vec::new();
        for name in &missing_seqids {
            let bytes = name.as_bytes();
            let len = u16_from_usize(bytes.len(), "missing seqid length")?;
            missing_seqid_table.extend_from_slice(&len.to_le_bytes());
            missing_seqid_table.extend_from_slice(bytes);
        }
        let missing_seqid_count = u32_from_usize(missing_seqids.len(), "missing seqid count")?;
        let missing_seqid_table_len =
            u32_from_usize(missing_seqid_table.len(), "missing seqid table length")?;

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

    /// Write one chrom block to disk and record its directory entry.
    /// The `ChromBlockBuilder` is consumed and its memory freed after this call.
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

        // write junction pool (zstd-compressed) next to tx_bytes
        let junction_pool_offset = tx_offset + tx_bytes as u64;
        let junction_pool_len = encode_compressed(&entry.junction_pool, &mut self.file)?;
        self.current_offset += junction_pool_len as u64;

        // then string pool (zstd-compressed)
        let string_pool_offset = junction_pool_offset + junction_pool_len as u64;
        let string_pool_len = encode_compressed(&entry.string_pool, &mut self.file)?;
        self.current_offset += string_pool_len as u64;

        // then splice site pool (zstd-compressed)
        let splice_site_pool_offset = string_pool_offset + string_pool_len as u64;
        let splice_site_pool_len = encode_compressed(&entry.splice_site_pool, &mut self.file)?;
        self.current_offset += splice_site_pool_len as u64;

        // generate a chromsome directory entry
        // and insert it into entries
        // each add_chrom will generate one entry on the fly
        let tx_count = u32_from_usize(entry.txs.len(), "chromosome transcript count")?;
        self.header.total_tx_n = self
            .header
            .total_tx_n
            .checked_add(tx_count)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "total transcript count exceeded u32",
                )
            })?;
        self.entries.push(ChromDirectoryEntry {
            chrom_id: entry.chrom_id,
            chrom_name_offset,
            chrom_name_len,
            global_tx_count: tx_count,
            global_tx_offset: u32_from_u64(tx_offset, "global tx offset")?,
            global_junction_pool_offset: u32_from_u64(
                junction_pool_offset,
                "global junction pool offset",
            )?,
            global_junction_count: junction_pool_len,
            global_string_pool_offset: u32_from_u64(
                string_pool_offset,
                "global string pool offset",
            )?,
            global_string_len: string_pool_len,
            global_splice_site_pool_offset: u32_from_u64(
                splice_site_pool_offset,
                "global splice site pool offset",
            )?,
            global_splice_site_pool_len: splice_site_pool_len,
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
