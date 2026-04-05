use std::collections::HashMap;

use crate::{
    core::tx_base_error::*,
    traits::{Encodable, PartialLoad},
};

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
