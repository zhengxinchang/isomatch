use crate::core::core_error::TxBaseError;
use crate::core::junction_pool::JunctionPool;
use crate::core::splice_site_pair::SpliceSitePair;
use crate::core::splice_site_pool::SpliceSitePool;
use crate::core::splice_site_span::SpliceSiteSpan;
use crate::core::string_pool::StringPool;
use crate::core::tx_base::TxBase;
use crate::core::tx_strand::ISOMSTRAND;
use crate::index::IndexStats;
use crate::index::fasta::FastaReader;
use crate::index::gtf::TxStructure;
use crate::index::index_error::IndexError;
use crate::traits::{Decodable, DiskSize, Encodable};
use crate::utils;
use std::io::{Error, ErrorKind, Read, Seek, SeekFrom, Write};

// Flag for index status
// bit 0, sequence from reference genome (0) or tx sequence (1)
// bit 1, gtf format, 0: plan text, 1: bgzipped
pub struct Flags {
    pub bits: u64,
}

impl Flags {
    // bit 0: 0 = gtf format, 1 = bgzipped
    const GTF_FORMAT_BIT: u64 = 1 << 0;
    // bit 1: index has ref hash
    const REF_HASH_BIT: u64 = 1 << 1;
    // bit 2: index has seq length
    const SEQ_HASH_BIT: u64 = 1 << 2;

    pub fn new() -> Self {
        Self { bits: 0 }
    }

    /// true = ref genome sequence hash is valid
    pub fn set_ref_hash(&mut self, has_ref_hash: bool) {
        if has_ref_hash {
            self.bits &= !Self::REF_HASH_BIT; // clear → ref
        } else {
            self.bits |= Self::REF_HASH_BIT; // set → tx seq
        }
    }

    /// Returns true if sequence is from reference genome
    pub fn get_ref_hash(&self) -> bool {
        self.bits & Self::REF_HASH_BIT == 0
    }

    /// true = seq hash is valid
    pub fn set_seq_hash(&mut self, has_seq_hash: bool) {
        if has_seq_hash {
            self.bits &= !Self::SEQ_HASH_BIT; //
        } else {
            self.bits |= Self::SEQ_HASH_BIT; //
        }
    }

    /// Returns true if sequence hash is valid
    /// true = valid, false = invalid
    pub fn get_seq_hash(&self) -> bool {
        self.bits & Self::SEQ_HASH_BIT == 0
    }

    /// true = gtf format, false = bgzipped
    pub fn set_gtf_format(&mut self, is_bgzipped: bool) {
        if is_bgzipped {
            self.bits |= Self::GTF_FORMAT_BIT; // set → bgzipped
        } else {
            self.bits &= !Self::GTF_FORMAT_BIT; // clear → plain text
        }
    }

    /// Returns true if GTF is bgzipped
    /// true = bgzipped, false = plain text
    pub fn get_gtf_format(&self) -> bool {
        self.bits & Self::GTF_FORMAT_BIT != 0
    }
}

pub struct IndexHeader {
    pub magic: [u8; 4],
    pub version: u8,
    pub flags: Flags,
    pub chrom_count: u32,
    pub gtf_size: u64,
    pub index_size: u64,
    pub md5: [u8; 16],
    /// Byte length of the chrom name table that immediately follows the directory.
    pub chrom_name_table_len: u32,
    pub reserved_to_4k: [u8; 4096 - 4 - 1 - 8 - 4 - 8 - 8 - 16 - 4], // 4043 bytes
}

impl IndexHeader {
    pub const CURRENT_VERSION: u8 = 1;

    pub fn new(
        chrom_count: u32,
        gtf_size: u64,
        index_size: u64,
        md5: [u8; 16],
        has_ref_hash: bool,
        has_seq_hash: bool,
        chrom_name_table_len: u32,
    ) -> Self {
        let mut flags = Flags::new();
        flags.set_ref_hash(has_ref_hash);
        flags.set_seq_hash(has_seq_hash);
        Self {
            magic: *b"ISOM",
            version: Self::CURRENT_VERSION,
            flags,
            chrom_count,
            gtf_size,
            index_size,
            md5,
            chrom_name_table_len,
            reserved_to_4k: [0u8; 4096 - 4 - 1 - 8 - 4 - 8 - 8 - 16 - 4],
        }
    }
}

impl DiskSize for IndexHeader {
    const DISK_SIZE: usize = 4096; // 固定 4KB 大小
}

impl Encodable for IndexHeader {
    type Error = std::io::Error;

    fn encode_to<W: Write>(&self, writer: &mut W) -> Result<usize, Self::Error> {
        let mut buf = Vec::new();
        buf.extend_from_slice(&self.magic);
        buf.push(self.version);
        buf.extend_from_slice(&self.flags.bits.to_le_bytes());
        buf.extend_from_slice(&self.chrom_count.to_le_bytes());
        buf.extend_from_slice(&self.gtf_size.to_le_bytes());
        buf.extend_from_slice(&self.index_size.to_le_bytes());
        buf.extend_from_slice(&self.md5);
        buf.extend_from_slice(&self.chrom_name_table_len.to_le_bytes());
        buf.extend_from_slice(&self.reserved_to_4k);
        writer.write_all(&buf)?;
        Ok(buf.len())
    }
}

impl Decodable for IndexHeader {
    type Error = std::io::Error;
    type Args = ();

    fn decode_from<R: Read + Seek>(reader: &mut R, _args: Self::Args) -> Result<Self, Self::Error> {
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if magic != *b"ISOM" {
            return Err(Error::new(
                ErrorKind::InvalidData,
                "Verify index file completeness failed. The index may be corrupted. Please rerun `isomatch index`",
            ));
        }

        let mut version_buf = [0u8; 1];
        reader.read_exact(&mut version_buf)?;
        let version = version_buf[0];
        if version != Self::CURRENT_VERSION {
            return Err(Error::new(
                ErrorKind::InvalidData,
                format!(
                    "unsupported index version: expected {}, got {}",
                    Self::CURRENT_VERSION,
                    version
                ),
            ));
        }

        let mut flags_buf = [0u8; 8];
        reader.read_exact(&mut flags_buf)?;
        let flags = u64::from_le_bytes(flags_buf);
        let flags = Flags { bits: flags };

        let mut chrom_count_buf = [0u8; 4];
        reader.read_exact(&mut chrom_count_buf)?;
        let chrom_count = u32::from_le_bytes(chrom_count_buf);

        let mut gtf_size_buf = [0u8; 8];
        reader.read_exact(&mut gtf_size_buf)?;
        let gtf_size = u64::from_le_bytes(gtf_size_buf);

        let mut index_size_buf = [0u8; 8];
        reader.read_exact(&mut index_size_buf)?;
        let index_size = u64::from_le_bytes(index_size_buf);

        let mut md5 = [0u8; 16];
        reader.read_exact(&mut md5)?;

        let mut chrom_name_table_len_buf = [0u8; 4];
        reader.read_exact(&mut chrom_name_table_len_buf)?;
        let chrom_name_table_len = u32::from_le_bytes(chrom_name_table_len_buf);

        // consume remaining reserved bytes to stay at 4 KB boundary
        let mut reserved_to_4k = [0u8; 4096 - 4 - 1 - 8 - 4 - 8 - 8 - 16 - 4];
        reader.read_exact(&mut reserved_to_4k)?;

        if index_size < Self::DISK_SIZE as u64 {
            return Err(Error::new(
                ErrorKind::InvalidData,
                format!(
                    "invalid index size in header: {} is smaller than header size {}",
                    index_size,
                    Self::DISK_SIZE
                ),
            ));
        }

        let next_pos = reader.stream_position()?;
        let actual_index_size = reader.seek(SeekFrom::End(0))?;
        reader.seek(SeekFrom::Start(next_pos))?;

        if actual_index_size != index_size {
            return Err(Error::new(
                ErrorKind::InvalidData,
                format!(
                    "index size mismatch: header says {}, actual file size is {}",
                    index_size, actual_index_size
                ),
            ));
        }

        Ok(Self {
            magic,
            version,
            flags,
            chrom_count,
            gtf_size,
            index_size,
            md5,
            chrom_name_table_len,
            reserved_to_4k,
        })
    }
}

pub struct ChromDirectoryEntry {
    pub chrom_id: u16,
    pub chrom_name_offset: u32,
    pub chrom_name_len: u32,
    pub global_tx_offset: u32,
    pub global_tx_count: u32,
    pub global_junction_pool_offset: u32,
    pub global_junction_count: u32,
    pub global_string_pool_offset: u32,
    pub global_string_len: u32,
    pub global_splice_site_pool_offset: u32,
    pub global_splice_site_pool_len: u32,
}

impl DiskSize for ChromDirectoryEntry {
    const DISK_SIZE: usize = 42;
}

impl Encodable for ChromDirectoryEntry {
    type Error = std::io::Error;

    fn encode_to<W: Write>(&self, writer: &mut W) -> Result<usize, Self::Error> {
        writer.write_all(&self.chrom_id.to_le_bytes())?;
        writer.write_all(&self.chrom_name_offset.to_le_bytes())?;
        writer.write_all(&self.chrom_name_len.to_le_bytes())?;
        writer.write_all(&self.global_tx_count.to_le_bytes())?;
        writer.write_all(&self.global_tx_offset.to_le_bytes())?;
        writer.write_all(&self.global_junction_pool_offset.to_le_bytes())?;
        writer.write_all(&self.global_junction_count.to_le_bytes())?;
        writer.write_all(&self.global_string_pool_offset.to_le_bytes())?;
        writer.write_all(&self.global_string_len.to_le_bytes())?;
        writer.write_all(&self.global_splice_site_pool_offset.to_le_bytes())?;
        writer.write_all(&self.global_splice_site_pool_len.to_le_bytes())?;
        Ok(Self::DISK_SIZE)
    }
}

impl Decodable for ChromDirectoryEntry {
    type Error = std::io::Error;
    type Args = ();

    fn decode_from<R: Read + Seek>(reader: &mut R, _args: Self::Args) -> Result<Self, Self::Error> {
        let mut buf = [0u8; ChromDirectoryEntry::DISK_SIZE];
        reader.read_exact(&mut buf)?;

        Ok(Self {
            chrom_id: u16::from_le_bytes(buf[0..2].try_into().unwrap()),
            chrom_name_offset: u32::from_le_bytes(buf[2..6].try_into().unwrap()),
            chrom_name_len: u32::from_le_bytes(buf[6..10].try_into().unwrap()),
            global_tx_count: u32::from_le_bytes(buf[10..14].try_into().unwrap()),
            global_tx_offset: u32::from_le_bytes(buf[14..18].try_into().unwrap()),
            global_junction_pool_offset: u32::from_le_bytes(buf[18..22].try_into().unwrap()),
            global_junction_count: u32::from_le_bytes(buf[22..26].try_into().unwrap()),
            global_string_pool_offset: u32::from_le_bytes(buf[26..30].try_into().unwrap()),
            global_string_len: u32::from_le_bytes(buf[30..34].try_into().unwrap()),
            global_splice_site_pool_offset: u32::from_le_bytes(buf[34..38].try_into().unwrap()),
            global_splice_site_pool_len: u32::from_le_bytes(buf[38..42].try_into().unwrap()),
        })
    }
}

/// Builder for constructing a single chrom's data block.
pub struct ChromBlockBuilder {
    pub chrom_id: u16,
    pub txs: Vec<TxBase>,
    pub junction_pool: JunctionPool,
    pub splice_site_pool: SpliceSitePool,
    pub string_pool: StringPool,
}

impl ChromBlockBuilder {
    pub fn init(chrom_id: u16) -> Self {
        Self {
            chrom_id,
            txs: Vec::new(),
            junction_pool: JunctionPool::new(),
            splice_site_pool: SpliceSitePool::new(),
            string_pool: StringPool::new(),
        }
    }

    pub fn add_tx(
        &mut self,
        gtf_tx: TxStructure,
        refr: &mut FastaReader,
        seqr: &mut Option<FastaReader>,
        stats: &mut IndexStats,
    ) -> Result<(), IndexError> {
        // let intron:Vec<u32> = gtf_tx.exons.iter().flat_map(|(e1,e2)|[*e1,*e2]).collect::<Vec<_>>();
        let intron: Vec<u32> = gtf_tx
            .exons
            .windows(2)
            .flat_map(|w| [w[0].1, w[1].0]) // [exon_n.end, exon_{n+1}.start, ...]
            .collect();

        let junction_span =
            self.junction_pool
                .add(&intron)
                .map_err(|e| IndexError::JunctionPoolAdd {
                    id: gtf_tx.tx_id.clone(),
                    reason: e.to_string(),
                })?;

        let tx_id_span =
            self.string_pool
                .add(gtf_tx.tx_id.as_str())
                .map_err(|e| IndexError::StringPoolAdd {
                    id: gtf_tx.tx_id.clone(),
                    reason: e.to_string(),
                })?;

        let gene_id_span = self.string_pool.add(gtf_tx.gene_id.as_str()).map_err(|e| {
            IndexError::StringPoolAdd {
                id: gtf_tx.gene_id.clone(),
                reason: e.to_string(),
            }
        })?;

        let (refhash, splice_site_pairs) = if gtf_tx.exons.len() == 1 {
            // if mono-exon, dont need to do this calculation
            (0, Vec::new())
        } else {
            let reference_seq = refr
                .fetch(
                    &gtf_tx.chrom,
                    gtf_tx.get_0based_start() as usize,
                    gtf_tx.end as usize,
                    true,
                )
                .map_err(|e| IndexError::FetchSeqFailed {
                    reason: format!(
                        "Can not fetch reference sequence for transcript {} on {}:{}-{}: {}",
                        gtf_tx.tx_id, gtf_tx.chrom, gtf_tx.start, gtf_tx.end, e
                    ),
                })?;

            let mut exon_offsets: Vec<(u32, u32)> = gtf_tx.get_0based_exon_relative_offset();

            let splice_sites_offsets: Vec<(usize, usize, usize, usize)> = exon_offsets
                .windows(2)
                .map(|e| {
                    (
                        e[0].1 as usize,
                        e[0].1 as usize + 2,
                        (e[1].0 as usize).saturating_sub(2),
                        e[1].0 as usize,
                    )
                })
                .collect();
            // println!("{:?}",gtf_tx.exons);
            // println!("{:?}",exon_offsets);
            // println!("{}",reference_seq.len());

            // shift the first and the last position into 3bp close to another side of the exon.
            // only for first and last exon
            // this ensure that transcripts that share same isoform structure (including the small exon shift).
            // will be correctly assinged same hash vlaue.

            let left_exon = exon_offsets
                .first_mut()
                .ok_or_else(|| IndexError::FetchSeqFailed {
                    reason: "Can not obtain the frst exon".to_string(),
                })?;
            if (left_exon.1 - left_exon.0) > 3 {
                left_exon.0 = left_exon.1 - 3;
            }

            let right_exon = exon_offsets
                .last_mut()
                .ok_or_else(|| IndexError::FetchSeqFailed {
                    reason: "Can not obain the last exon".to_string(),
                })?;

            if (right_exon.1 - right_exon.0) > 3 {
                right_exon.1 = right_exon.0 + 3;
            }

            let mut tx_sequence = Vec::new();

            for region in exon_offsets.into_iter() {
                let bases = &reference_seq[region.0 as usize..region.1 as usize];
                tx_sequence.extend_from_slice(&bases);
            }

            let mut splice_site_pairs = Vec::new();
            for (lstart, lend, rstart, rend) in splice_sites_offsets.into_iter() {
                let left = reference_seq[lstart..lend].to_vec();
                let right = reference_seq[rstart..rend].to_vec();
                splice_site_pairs.push((left, right));
            }

            // refhash
            (utils::hash_u8_vec(&tx_sequence), splice_site_pairs)
        };

        // seqhash
        let seqhash = match seqr {
            Some(reader) => {
                if gtf_tx.exons.len() == 1 {
                    // if mono-exon, dont need to do this calculation
                    0
                } else {
                    let sequence = reader.fetch_all(&gtf_tx.tx_id, false)?;

                    let tx_seq_len: usize = gtf_tx
                        .exons
                        .iter()
                        .map(|(s, e)| (*e - *s + 1) as usize)
                        .sum();

                    if sequence.len() != tx_seq_len {
                        dbg!(&gtf_tx);
                        return Err(IndexError::FetchSeqFailed {
                            reason: format!("Actual sequence length ({}) is not equal to GTF derived sequence length ({}). Affected transcript {}",sequence.len(),tx_seq_len,gtf_tx.tx_id).to_string()
                        });
                    }

                    let first_exon = gtf_tx.exons.first().ok_or(IndexError::FetchSeqFailed {
                        reason: "No exons found".to_string(),
                    })?;

                    let last_exon = gtf_tx.exons.last().ok_or(IndexError::FetchSeqFailed {
                        reason: "No exons found".to_string(),
                    })?;

                    let first_exon_len = (first_exon.1 - first_exon.0 + 1) as usize;
                    let last_exon_len = (last_exon.1 - last_exon.0 + 1) as usize;

                    // exons are sorted by genomic position.
                    // For minus strand the RNA 5'-terminal exon is genomically last,
                    // and the RNA 3'-terminal exon is genomically first.
                    let (tss_exon_len, tes_exon_len) = if gtf_tx.strand == ISOMSTRAND::Minus {
                        (last_exon_len, first_exon_len)
                    } else {
                        (first_exon_len, last_exon_len)
                    };

                    // trim from the outer (TSS) side of the 5'-terminal exon
                    let left_trim = if tss_exon_len > 3 {
                        tss_exon_len - 3
                    } else {
                        0
                    };

                    // trim from the outer (TES) side of the 3'-terminal exon
                    let right_trim = if tes_exon_len > 3 {
                        sequence.len() - (tes_exon_len - 3)
                    } else {
                        sequence.len()
                    };

                    let sliced_seq = &sequence[left_trim..right_trim];

                    utils::hash_u8_slice(&sliced_seq)
                }
            }
            None => 0u128,
        };

        let canonical_junction_count = splice_site_pairs
            .iter()
            .map(|(left, right)| SpliceSitePair::pack(left, right, gtf_tx.strand))
            .try_fold(0usize, |count, pair| {
                let pair: SpliceSitePair = pair?;
                Ok::<usize, TxBaseError>(if pair.is_canonical() {
                    count + 1
                } else {
                    count
                })
            })
            .map_err(|e| IndexError::AddGTFTx {
                id: gtf_tx.tx_id.clone(),
                reason: e.to_string(),
            })?;

        // build splice site pairs from intron boundaries
        // let splice_site_pairs: Vec<(&str, &str)> = Vec::new(); // TODO: extract donor/acceptor dinucleotides from reference
        let splice_site_span = if splice_site_pairs.is_empty() {
            SpliceSiteSpan {
                offset: 0,
                count: 0,
            }
        } else {
            self.splice_site_pool
                .add_pairs(&splice_site_pairs, gtf_tx.strand)
                .map_err(|e| IndexError::AddGTFTx {
                    id: gtf_tx.tx_id.clone(),
                    reason: e.to_string(),
                })?
        };

        let tx_base = TxBase::new(
            gtf_tx.idx,
            self.chrom_id,
            gtf_tx.start,
            gtf_tx.end,
            gtf_tx.strand,
            seqhash,
            refhash,
            gtf_tx.exons.len() as u16,
            splice_site_span,
            junction_span,
            tx_id_span,
            gene_id_span,
        )
        .map_err(|e| IndexError::AddGTFTx {
            id: gtf_tx.tx_id.clone(),
            reason: e.to_string(),
        })?;

        stats.observe_tx(
            gtf_tx.strand,
            gtf_tx.exons.len(),
            canonical_junction_count,
            &gtf_tx.gene_id,
        );
        self.txs.push(tx_base);

        Ok(())
    }

    pub fn finalize(&mut self) {
        self.txs.sort_unstable();
    }
}
