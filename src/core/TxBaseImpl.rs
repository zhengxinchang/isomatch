use crate::core::TxBase::TxBaseFlags;
use crate::core::TxBase::{JunctionSpan, StringSpan, TxBase, TxBaseTrait};
use crate::core::TxBaseError::TxBaseError;
use crate::core::TxBoundary::TxBoundary;
use crate::traits::{DiskSize, Encodable, PartialLoad};

impl TxBaseTrait for TxBase {
    fn tx_idx(&self) -> u32 {
        self.tx_idx
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

    fn gtf_offset(&self) -> u64 {
        self._gtf_offset
    }

    fn gtf_len(&self) -> u32 {
        self._gtf_len
    }

    fn n_exons(&self) -> u16 {
        self.n_exons
    }

    fn junctions(&self) -> JunctionSpan {
        self.junctions
    }

    fn transcript_span(&self) -> StringSpan {
        self.tx_id_span
    }

    fn gene_span(&self) -> StringSpan {
        self.gene_id_span
    }
}

impl DiskSize for TxBase {
    const DISK_SIZE: usize = 84;
}

impl Encodable for TxBase {
    type Error = TxBaseError;

    /// no need to do the TxBoundary encoding here since TxBase already stores start, end, strand separately for easy access
    fn encode_to<W: std::io::Write>(&self, writer: &mut W) -> Result<usize, Self::Error> {
        writer
            .write_all(&self.tx_idx.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.chrom_id.to_le_bytes())
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
            .write_all(&self._gtf_offset.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self._gtf_len.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.n_exons.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.junctions.offset.to_le_bytes())
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        writer
            .write_all(&self.junctions.count.to_le_bytes())
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
    type Args = (); // TxBase is self-contained, no extra context needed

    fn load_range<R: std::io::Read + std::io::Seek>(
        reader: &mut R,
        offset: u64,
        _len: usize, // always DISK_SIZE for fixed-size TxBase, ignored
        _args: Self::Args,
    ) -> Result<Self, Self::Error> {
        reader
            .seek(std::io::SeekFrom::Start(offset))
            .map_err(|e| TxBaseError::Io(e.to_string()))?;

        let mut buf = [0u8; TxBase::DISK_SIZE];
        reader
            .read_exact(&mut buf)
            .map_err(|e| TxBaseError::Io(e.to_string()))?;

        let tx_id = u32::from_le_bytes(buf[0..4].try_into().unwrap());
        let chrom_id = u16::from_le_bytes(buf[4..6].try_into().unwrap());
        let start = u32::from_le_bytes(buf[6..10].try_into().unwrap());
        let end = u32::from_le_bytes(buf[10..14].try_into().unwrap());
        let flags = TxBaseFlags(u16::from_le_bytes(buf[14..16].try_into().unwrap()));
        let seq_hash = u128::from_le_bytes(buf[16..32].try_into().unwrap());
        let ref_hash = u128::from_le_bytes(buf[32..48].try_into().unwrap());
        let gtf_offset = u64::from_le_bytes(buf[48..56].try_into().unwrap());
        let gtf_len = u32::from_le_bytes(buf[56..60].try_into().unwrap());
        let n_exons = u16::from_le_bytes(buf[60..62].try_into().unwrap());
        let junctions_offset = u32::from_le_bytes(buf[62..66].try_into().unwrap());
        let junctions_count = u16::from_le_bytes(buf[66..68].try_into().unwrap());
        let transcript_span_offset = u32::from_le_bytes(buf[68..72].try_into().unwrap());
        let transcript_span_byte_len = u32::from_le_bytes(buf[72..76].try_into().unwrap());
        let gene_span_offset = u32::from_le_bytes(buf[76..80].try_into().unwrap());
        let gene_span_byte_len = u32::from_le_bytes(buf[80..84].try_into().unwrap());

        Ok(Self {
            tx_idx: tx_id,
            boundary: TxBoundary::new(start, end, flags.strand()),
            chrom_id,
            start,
            end,
            flags,
            seq_hash,
            ref_hash,
            _gtf_offset: gtf_offset,
            _gtf_len: gtf_len,
            n_exons,
            junctions: JunctionSpan {
                offset: junctions_offset,
                count: junctions_count,
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
            tx_idx: 7,
            boundary: TxBoundary::new(101, 250, 1),
            chrom_id: 3,
            start: 101,
            end: 250,
            flags: TxBaseFlags::new(1).unwrap(),
            seq_hash: 11,
            ref_hash: 22,
            _gtf_offset: 1234,
            _gtf_len: 56,
            n_exons: 2,
            junctions: JunctionSpan {
                offset: 9,
                count: 2,
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

        let decoded = TxBase::load_range(&mut Cursor::new(buf), 0, TxBase::DISK_SIZE, ()).unwrap();

        assert_eq!(decoded.tx_idx, tx.tx_idx);
        assert_eq!(decoded.seq_hash, tx.seq_hash);
        assert_eq!(decoded.ref_hash, tx.ref_hash);
        assert_eq!(decoded.gtf_offset(), tx.gtf_offset());
        assert_eq!(decoded.tx_id_span, tx.tx_id_span);
        assert_eq!(decoded.gene_id_span, tx.gene_id_span);
    }
}
