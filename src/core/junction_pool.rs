use crate::{
    core::tx_base_error::TxBaseError,
    traits::{Encodable, PartialLoad},
};

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
        if junctions.len() % 2 != 0 {
            return Err(TxBaseError::InvalidEncoding {
                msg: format!(
                    "junction coordinate count must be even, got {}",
                    junctions.len()
                ),
            });
        }

        if junctions.chunks_exact(2).any(|pair| pair[0] >= pair[1])
            || junctions.windows(2).any(|window| window[0] > window[1])
        {
            return Err(TxBaseError::JunctionsNotStrictlyIncreasing);
        }
        Ok(())
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
