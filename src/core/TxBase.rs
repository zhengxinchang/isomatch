use std::collections::HashMap;

use crate::{
    core::{TxBaseError::TxBaseError, TxBoundary::TxBoundary},
    traits::{Encodable, PartialLoad},
};

pub trait TxBaseTrait {
    fn tx_idx(&self) -> u32;
    fn tx_boundary(&self) -> TxBoundary {
        TxBoundary::new(self.start(), self.end(), self.strand())
    }
    fn chrom_id(&self) -> u16;
    fn start(&self) -> u32;
    fn end(&self) -> u32;
    fn flags(&self) -> TxBaseFlags;
    fn seq_hash(&self) -> u128;
    fn ref_hash(&self) -> u128;
    fn gtf_offset(&self) -> u64;
    fn gtf_len(&self) -> u32;
    fn n_exons(&self) -> u16;
    fn junctions(&self) -> JunctionSpan;
    fn transcript_span(&self) -> StringSpan;
    fn gene_span(&self) -> StringSpan;
    fn strand(&self) -> u8 {
        self.flags().strand()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TxBase {
    pub tx_idx: u32,
    pub boundary: TxBoundary,
    pub chrom_id: u16,
    pub start: u32,
    pub end: u32,
    pub flags: TxBaseFlags,
    pub seq_hash: u128,
    pub ref_hash: u128,
    pub _gtf_offset: u64, // byte offset of the GTF record in the original GTF file
    pub _gtf_len: u32,    // byte length of the GTF record in the original GTF file
    pub n_exons: u16,
    pub junctions: JunctionSpan,
    /// Direct reference into the on-disk string section for GTF `transcript_id`.
    pub tx_id_span: StringSpan,
    /// Direct reference into the on-disk string section for GTF `gene_id`.
    pub gene_id_span: StringSpan,
}

impl TxBase {
    pub fn new(
        tx_idx: u32, // record index in the GTF file
        chrom_id: u16,
        start: u32,
        end: u32,
        strand: u8,
        seq_hash: u128,
        ref_hash: u128,
        n_exons: u16,
        junctions: JunctionSpan,
        transcript_span: StringSpan,
        gene_span: StringSpan,
    ) -> Result<Self, TxBaseError> {
        if start > end {
            return Err(TxBaseError::InvalidBounds { start, end });
        }
        if n_exons == 0 {
            return Err(TxBaseError::InvalidExonCount { n_exons });
        }

        Ok(Self {
            tx_idx: tx_idx,
            boundary: TxBoundary::new(start, end, strand),
            chrom_id,
            start,
            end,
            flags: TxBaseFlags::new(strand)?,
            seq_hash,
            ref_hash,
            _gtf_offset: 0,
            _gtf_len: 0,
            n_exons,
            junctions,
            tx_id_span: transcript_span,
            gene_id_span: gene_span,
        })
    }

    pub fn strand(&self) -> u8 {
        self.flags.strand()
    }

    pub fn sort_key(&self) -> (u16, u32, u32, u8) {
        (self.chrom_id, self.start, self.end, self.strand())
    }

    pub fn junction_slice<'a>(&self, pool: &'a JunctionPool) -> Result<&'a [u32], TxBaseError> {
        pool.get(self.junctions)
    }
}

/// Flags for TxBase.
/// bit 0: strand (0 for +, 1 for -)
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
#[repr(transparent)]
pub struct TxBaseFlags(pub u16);

impl TxBaseFlags {
    const NEG_STRAND_BIT: u16 = 1;

    pub fn new(strand: u8) -> Result<Self, TxBaseError> {
        match strand {
            0 => Ok(Self(0)),
            1 => Ok(Self(Self::NEG_STRAND_BIT)),
            _ => Err(TxBaseError::InvalidStrand { strand }),
        }
    }

    pub fn strand(self) -> u8 {
        if self.0 & Self::NEG_STRAND_BIT == Self::NEG_STRAND_BIT {
            1
        } else {
            0
        }
    }

    pub fn bits(self) -> u16 {
        self.0
    }

    /// Set bit at position `bit` (0-indexed), returns new Flags.
    pub fn set_bit(self, bit: u16) -> Self {
        Self(self.0 | (1u16 << bit))
    }

    /// Clear bit at position `bit` (0-indexed), returns new Flags.
    pub fn clear_bit(self, bit: u16) -> Self {
        Self(self.0 & !(1u16 << bit))
    }

    /// Get the value of bit at position `bit` (0-indexed).
    pub fn get_bit(self, bit: u16) -> bool {
        (self.0 >> bit) & 1 == 1
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct JunctionSpan {
    pub offset: u32,
    pub count: u16,
}

impl JunctionSpan {
    pub fn is_empty(self) -> bool {
        self.count == 0
    }

    pub fn end_offset(self) -> u32 {
        self.offset + u32::from(self.count)
    }
}

#[derive(Debug)]
pub struct JunctionPool {
    pub junctions: Vec<u32>,
    // pub disk_handle: Option<BufReader<std::fs::File>>,
}

impl JunctionPool {
    const CANONICAL_BIT: u32 = 1;
    const COORD_MASK: u32 = !Self::CANONICAL_BIT;

    pub fn new() -> Self {
        Self {
            junctions: Vec::new(),
            // disk_handle: None,
        }
    }

    pub fn with_capacity(capacity: usize) -> Result<Self, TxBaseError> {
        Ok(Self {
            junctions: Vec::with_capacity(capacity),
            // disk_handle: None,
        })
    }

    pub fn add(&mut self, junctions: &Vec<u32>) -> Result<JunctionSpan, TxBaseError> {
        Self::validate_junctions(junctions)?;

        let offset = u32::try_from(self.junctions.len()).map_err(|_| TxBaseError::PoolTooLarge)?;
        let count = u16::try_from(junctions.len()).map_err(|_| TxBaseError::TooManyJunctions {
            count: junctions.len(),
        })?;

        self.junctions.extend_from_slice(junctions);

        Ok(JunctionSpan { offset, count })
    }

    pub fn get(&self, span: JunctionSpan) -> Result<&[u32], TxBaseError> {
        let start: usize = usize::try_from(span.offset).map_err(|_| TxBaseError::InvalidSpan {
            offset: span.offset,
            count: span.count,
            pool_len: self.junctions.len(),
        })?;
        let end: usize = start + usize::from(span.count);

        self.junctions
            .get(start..end)
            .ok_or(TxBaseError::InvalidSpan {
                offset: span.offset,
                count: span.count,
                pool_len: self.junctions.len(),
            })
    }

    pub fn len(&self) -> usize {
        self.junctions.len()
    }

    pub fn is_empty(&self) -> bool {
        self.junctions.is_empty()
    }

    fn validate_junctions(junctions: &[u32]) -> Result<(), TxBaseError> {
        if junctions.windows(2).any(|window| window[0] >= window[1]) {
            return Err(TxBaseError::JunctionsNotStrictlyIncreasing);
        }
        Ok(())
    }

    pub fn encode_site(coord: u32, is_canonical: bool) -> u32 {
        (coord << 1) | (is_canonical as u32)
    }

    pub fn decode_coord(raw: u32) -> u32 {
        raw >> 1
    }

    pub fn decode_canonical(raw: u32) -> bool {
        (raw & Self::CANONICAL_BIT) == 1
    }
}

impl Encodable for JunctionPool {
    type Error = TxBaseError;
    fn encode_to<W: std::io::Write>(&self, writer: &mut W) -> Result<usize, Self::Error> {
        let mut written = 0;
        for &j in &self.junctions {
            writer
                .write_all(&j.to_le_bytes())
                .map_err(|e| TxBaseError::Io(e.to_string()))?;
            written += 4;
        }
        Ok(written)
    }
}

impl PartialLoad for JunctionPool {
    type Error = TxBaseError;
    type Args = u16; // chrom_id
    fn load_range<R: std::io::Read + std::io::Seek>(
        reader: &mut R,
        offset: u64,
        len: usize,
        _chrom_id: Self::Args,
    ) -> Result<Self, Self::Error> {
        let mut buf = vec![0; len];
        reader
            .seek(std::io::SeekFrom::Start(offset))
            .map_err(TxBaseError::io)?;
        reader
            .read_exact(&mut buf)
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        if buf.len() % 4 != 0 {
            return Err(TxBaseError::InvalidEncoding {
                msg: format!("junction data length {} is not a multiple of 4", buf.len()),
            });
        }
        let junctions = buf
            .chunks_exact(4)
            .map(|c| u32::from_le_bytes(c.try_into().unwrap()))
            .collect();
        Ok(Self { junctions })
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct StringSpan {
    /// Byte offset within the string-data section of the index file.
    pub offset: u32,
    /// Length in bytes of the UTF-8 encoded string.
    pub byte_len: u32,
}

impl StringSpan {
    pub const EMPTY: Self = Self {
        offset: 0,
        byte_len: 0,
    };

    pub fn is_empty(self) -> bool {
        self.byte_len == 0
    }
}

#[derive(Debug, Default)]
pub struct StringPool {
    strings: Vec<u8>,
    index: HashMap<String, StringSpan>,
}

impl StringPool {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add(&mut self, s: &str) -> Result<StringSpan, TxBaseError> {
        if self.index.len() >= u32::MAX as usize {
            return Err(TxBaseError::StringPoolTooLarge);
        }
        if self.index.contains_key(s) {
            return Ok(*self.index.get(s).unwrap());
        } else {
            // No this string found in the pool, add it to the pool
            let offset = self.strings.len() as u32;
            let byte_len = s.len() as u32;
            self.strings.extend_from_slice(s.as_bytes());
            let span = StringSpan { offset, byte_len };
            self.index.insert(s.to_string(), span);
            Ok(span)
        }
    }

    pub fn get(&self, span: StringSpan) -> Result<&str, TxBaseError> {
        let offset = usize::try_from(span.offset)
            .map_err(|_| TxBaseError::InvalidInternId { id: span.offset })?;
        let end = offset
            + usize::try_from(span.byte_len)
                .map_err(|_| TxBaseError::InvalidInternId { id: span.offset })?;
        let bytes = self
            .strings
            .get(offset..end)
            .ok_or(TxBaseError::InvalidInternId { id: span.offset })?;
        std::str::from_utf8(bytes).map_err(|_| TxBaseError::InvalidInternId { id: span.offset })
    }
}

impl Encodable for StringPool {
    type Error = TxBaseError;
    fn encode_to<W: std::io::Write>(&self, writer: &mut W) -> Result<usize, Self::Error> {
        writer
            .write_all(&self.strings)
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        Ok(self.strings.len())
    }
}

impl PartialLoad for StringPool {
    type Error = TxBaseError;
    type Args = ();
    fn load_range<R: std::io::Read + std::io::Seek>(
        reader: &mut R,
        offset: u64,
        len: usize,
        _args: Self::Args,
    ) -> Result<Self, Self::Error> {
        let mut buf = vec![0; len];
        reader
            .seek(std::io::SeekFrom::Start(offset))
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        reader
            .read_exact(&mut buf)
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        Ok(Self {
            strings: buf,
            index: HashMap::new(), // Index is not needed for partial loads
        })
    }
}
