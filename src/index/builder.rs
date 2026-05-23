use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};

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
    let len = compressed.len() as u32;
    writer.write_all(&compressed)?;
    Ok(len)
}

use crate::core::tx_base::TxBase;
use crate::index::format::{ChromBlockBuilder, ChromDirectoryEntry, IndexHeader};
use crate::traits::{DiskSize, Encodable};

/// Index 写入器。
///
/// # 写入流程
/// 1. `new()` — 写入 4KB 占位 Header + N×34B 占位 Directory + Chrom Name Table + Missing SeqID Table（一次写定）
/// 2. `add_chrom()` — 逐 chrom 顺序写入数据，写完立即可 drop，记录真实 offset 到 `entries`
/// 3. `finalize()` — seek 回文件头，回填真实 Header 与 Directory
pub struct IndexBuilder {
    header: IndexHeader,
    entries: Vec<ChromDirectoryEntry>,
    /// Pre-computed (offset_in_table, len) for each chrom, indexed by chrom_id - 1.
    chrom_name_offsets: Vec<(u32, u32)>,
    current_offset: u64,
    file: BufWriter<File>,
}

impl IndexBuilder {
    /// `chrom_names` must be in the same order as chroms will appear in the GTF
    /// (i.e. the order returned by `profile_gtf`).
    pub fn new(
        file: File,
        chrom_names: Vec<String>,
        gtf_size: u64,
        md5: [u8; 16],
        has_ref_hash: bool,
        has_seq_hash: bool,
        missing_seqids: Vec<String>,
    ) -> std::io::Result<Self> {
        let mut file = BufWriter::new(file);
        let chrom_count = chrom_names.len() as u32;

        // Build the chrom name table bytes and pre-compute per-chrom offsets.
        let mut name_table: Vec<u8> = Vec::new();
        let mut chrom_name_offsets: Vec<(u32, u32)> = Vec::with_capacity(chrom_names.len());
        for name in &chrom_names {
            let offset = name_table.len() as u32;
            let len = name.len() as u32;
            name_table.extend_from_slice(name.as_bytes());
            chrom_name_offsets.push((offset, len));
        }
        let chrom_name_table_len = name_table.len() as u32;

        // Build the missing seqid table: each entry is u16 len + utf-8 bytes.
        let mut missing_seqid_table: Vec<u8> = Vec::new();
        for name in &missing_seqids {
            let bytes = name.as_bytes();
            missing_seqid_table.extend_from_slice(&(bytes.len() as u16).to_le_bytes());
            missing_seqid_table.extend_from_slice(bytes);
        }
        let missing_seqid_count = missing_seqids.len() as u32;
        let missing_seqid_table_len = missing_seqid_table.len() as u32;

        let header = IndexHeader::new(
            chrom_count,
            gtf_size,
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
        self.entries.push(ChromDirectoryEntry {
            chrom_id: entry.chrom_id,
            chrom_name_offset,
            chrom_name_len,
            global_tx_count: entry.txs.len() as u32,
            global_tx_offset: tx_offset as u32,
            global_junction_pool_offset: junction_pool_offset as u32,
            global_junction_count: junction_pool_len,
            global_string_pool_offset: string_pool_offset as u32,
            global_string_len: string_pool_len,
            global_splice_site_pool_offset: splice_site_pool_offset as u32,
            global_splice_site_pool_len: splice_site_pool_len,
        });

        Ok(())
    }

    /// Seek back and write the real header and directory.
    pub fn finalize(mut self) -> std::io::Result<()> {
        self.header.index_size = self.current_offset;
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
