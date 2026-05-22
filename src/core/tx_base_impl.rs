use crate::core::core_error::TxBaseError;
use crate::core::junction_pool::{JunctionPool, JunctionSpan};
use crate::core::splice_site_pair::SpliceSitePair;
use crate::core::splice_site_pool::SpliceSitePool;
use crate::core::splice_site_span::SpliceSiteSpan;
use crate::core::string_pool::{StringPool, StringSpan};
use crate::core::tx_base::{TxBase, TxBaseTrait};
use crate::core::tx_base_flag::TxBaseFlags;
use crate::core::tx_boundary::TxBoundary;
use crate::traits::{DiskSize, Encodable, PartialLoad};

#[derive(Debug, Clone, Copy)]
pub struct TxBaseLoadArgs {
    pub chrom_id: u16,
}

impl TxBaseTrait for TxBase {
    // 标注一下这个不应该被使用
    fn tx_idx(&self) -> u32 {
        self._tx_idx
    }
    fn tx_boundary(&self) -> TxBoundary {
        self.boundary
    }
    fn chrom_id(&self) -> u16 {
        self.chrom_id
    }

    fn start(&self) -> u32 {
        self.start
    }

    fn end(&self) -> u32 {
        self.end
    }

    fn flags(&self) -> TxBaseFlags {
        self.flags
    }

    fn seq_hash(&self) -> u128 {
        self.seq_hash
    }

    fn ref_hash(&self) -> u128 {
        self.ref_hash
    }

    // fn gtf_offset(&self) -> u64 {
    //     0
    // }

    // fn gtf_len(&self) -> u32 {
    //     0
    // }

    fn n_exons(&self) -> u16 {
        self.n_exons
    }

    fn junctions(&self, junction_pool: &JunctionPool, string_pool: &StringPool) -> Vec<(u32, u32)> {
        let raw = junction_pool.get(self.junctions_span).unwrap_or_else(|e| {
            panic!(
                "failed to resolve junction span for transcript {}: {}",
                self.source_tx_id(string_pool),
                e
            )
        });

        let chunks = raw.chunks_exact(2);
        let remainder = chunks.remainder();
        assert!(
            remainder.is_empty(),
            "junction coordinate count must be even for transcript {}, got {}",
            self.source_tx_id(string_pool),
            raw.len()
        );

        chunks.map(|pair| (pair[0], pair[1])).collect()
    }

    fn splice_sites(
        &self,
        splice_sites_pool: &SpliceSitePool,
        string_pool: &StringPool,
    ) -> Vec<SpliceSitePair> {
        splice_sites_pool
            .get_pair(self.splice_sites_span)
            .unwrap_or_else(|e| {
                panic!(
                    "failed to resolve splice site span for transcript {}: {}",
                    self.source_tx_id(string_pool),
                    e
                )
            })
            .to_vec()
    }

    fn source_tx_id(&self, string_pool: &StringPool) -> String {
        string_pool
            .get(self.tx_id_span)
            .unwrap_or_else(|e| {
                panic!(
                    "failed to resolve transcript_id span for transcript {}: {}",
                    self.source_tx_id(string_pool),
                    e
                )
            })
            .to_owned()
    }

    fn source_gene_id(&self, string_pool: &StringPool) -> String {
        string_pool
            .get(self.gene_id_span)
            .unwrap_or_else(|e| {
                panic!(
                    "failed to resolve gene_id span for transcript {}: {}",
                    self.source_tx_id(string_pool),
                    e
                )
            })
            .to_owned()
    }
}

impl DiskSize for TxBase {
    // Current layout drops on-disk chrom_id and the unused gtf_offset/gtf_len fields.
    const DISK_SIZE: usize = 76;
}

impl Encodable for TxBase {
    type Error = TxBaseError;

    /// no need to do the TxBoundary encoding here since TxBase already stores start, end, strand separately for easy access
    fn encode_to<W: std::io::Write>(&self, writer: &mut W) -> Result<usize, Self::Error> {
        writer
            .write_all(&self._tx_idx.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.start.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.end.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.flags.bits().to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.seq_hash.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.ref_hash.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.n_exons.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.junctions_span.offset.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.junctions_span.count.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.splice_sites_span.offset.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.splice_sites_span.count.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.tx_id_span.offset.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.tx_id_span.byte_len.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.gene_id_span.offset.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.gene_id_span.byte_len.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        Ok(Self::DISK_SIZE)
    }
}

impl PartialLoad for TxBase {
    type Error = TxBaseError;
    type Args = TxBaseLoadArgs;

    fn load_range<R: std::io::Read + std::io::Seek>(
        reader: &mut R,
        offset: u64,
        _len: usize, // always DISK_SIZE for fixed-size TxBase, ignored
        args: Self::Args,
    ) -> Result<Self, Self::Error> {
        reader
            .seek(std::io::SeekFrom::Start(offset))
            .map_err(|e| TxBaseError::Io(e.to_string()))?;

        let mut buf = [0u8; TxBase::DISK_SIZE];
        reader
            .read_exact(&mut buf)
            .map_err(|e| TxBaseError::Io(e.to_string()))?;

        let tx_id = u32::from_le_bytes(buf[0..4].try_into().unwrap());
        let start = u32::from_le_bytes(buf[4..8].try_into().unwrap());
        let end = u32::from_le_bytes(buf[8..12].try_into().unwrap());
        let flags = TxBaseFlags(u16::from_le_bytes(buf[12..14].try_into().unwrap()));
        let seq_hash = u128::from_le_bytes(buf[14..30].try_into().unwrap());
        let ref_hash = u128::from_le_bytes(buf[30..46].try_into().unwrap());
        let n_exons = u16::from_le_bytes(buf[46..48].try_into().unwrap());
        let junctions_offset = u32::from_le_bytes(buf[48..52].try_into().unwrap());
        let junctions_count = u16::from_le_bytes(buf[52..54].try_into().unwrap());
        let splice_sites_offset = u32::from_le_bytes(buf[54..58].try_into().unwrap());
        let splice_sites_count = u16::from_le_bytes(buf[58..60].try_into().unwrap());
        let transcript_span_offset = u32::from_le_bytes(buf[60..64].try_into().unwrap());
        let transcript_span_byte_len = u32::from_le_bytes(buf[64..68].try_into().unwrap());
        let gene_span_offset = u32::from_le_bytes(buf[68..72].try_into().unwrap());
        let gene_span_byte_len = u32::from_le_bytes(buf[72..76].try_into().unwrap());

        Ok(Self {
            _tx_idx: tx_id,
            boundary: TxBoundary::new(start, end, flags.get_strand()),
            chrom_id: args.chrom_id,
            start,
            end,
            flags,
            seq_hash,
            ref_hash,
            n_exons,
            junctions_span: JunctionSpan {
                offset: junctions_offset,
                count: junctions_count,
            },
            splice_sites_span: SpliceSiteSpan {
                offset: splice_sites_offset,
                count: splice_sites_count,
            },
            tx_id_span: StringSpan {
                offset: transcript_span_offset,
                byte_len: transcript_span_byte_len,
            },
            gene_id_span: StringSpan {
                offset: gene_span_offset,
                byte_len: gene_span_byte_len,
            },
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn txbase_roundtrip_preserves_ref_hash() {
        let tx = TxBase {
            _tx_idx: 7,
            boundary: TxBoundary::new(101, 250, crate::core::tx_strand::ISOMSTRAND::Minus),
            chrom_id: 3,
            start: 101,
            end: 250,
            flags: TxBaseFlags::new(crate::core::tx_strand::ISOMSTRAND::Minus, true).unwrap(),
            seq_hash: 11,
            ref_hash: 22,
            n_exons: 2,
            junctions_span: JunctionSpan {
                offset: 9,
                count: 2,
            },
            splice_sites_span: SpliceSiteSpan {
                offset: 0,
                count: 0,
            },
            tx_id_span: StringSpan {
                offset: 100,
                byte_len: 12,
            },
            gene_id_span: StringSpan {
                offset: 200,
                byte_len: 8,
            },
        };

        let mut buf = Vec::new();
        tx.encode_to(&mut buf).unwrap();
        assert_eq!(buf.len(), TxBase::DISK_SIZE);

        let decoded = TxBase::load_range(
            &mut Cursor::new(buf),
            0,
            TxBase::DISK_SIZE,
            TxBaseLoadArgs {
                chrom_id: tx.chrom_id,
            },
        )
        .unwrap();

        assert_eq!(decoded._tx_idx, tx._tx_idx);
        assert_eq!(decoded.chrom_id, tx.chrom_id);
        assert_eq!(decoded.seq_hash, tx.seq_hash);
        assert_eq!(decoded.ref_hash, tx.ref_hash);
        // assert_eq!(decoded.gtf_offset(), 0);
        assert_eq!(decoded.tx_id_span, tx.tx_id_span);
        assert_eq!(decoded.gene_id_span, tx.gene_id_span);
    }
}
