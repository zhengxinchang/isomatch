use std::io::{Read, Seek, Write};

/// Types with a fixed, known byte size on disk.
///
/// Implement this for every fixed-size record (`TxBase`, `JunctionSpan`,
/// `StringSpan`, `Flags`, etc.) so that writers can compute offsets and
/// readers can seek directly to record N via `base + N * DISK_SIZE`.
pub trait DiskSize {
    const DISK_SIZE: usize;
}

/// Serialize `self` into a writer.
///
/// For fixed-size types this writes exactly [`DiskSize::DISK_SIZE`] bytes.
/// For variable-length types (pools) it writes the full payload.
pub trait Encodable {
    type Error;
    fn encode_to<W: Write>(&self, writer: &mut W) -> Result<usize, Self::Error>;

    /// Convenience: encode into a new `Vec<u8>`.
    fn encode(&self) -> Result<Vec<u8>, Self::Error> {
        let mut buf = Vec::new();
        self.encode_to(&mut buf)?;
        Ok(buf)
    }
}

/// Reconstruct a value by reading from a byte source.
///
/// Uses an associated `Args` type for any context needed at decode time
/// (e.g. `chrom_id` for `JunctionPool`).  Use `()` when no extra
/// context is required.
pub trait Decodable: Sized {
    type Error;
    type Args;
    fn decode_from<R: Read + Seek>(reader: &mut R, args: Self::Args) -> Result<Self, Self::Error>;
}

/// Read a portion of a section from a seekable source.
///
/// Useful for on-demand loading of junction pools or string data
/// without reading the entire section into memory.
pub trait PartialLoad: Sized {
    type Error;
    type Args;
    fn load_range<R: Read + Seek>(
        reader: &mut R,
        offset: u64,
        len: usize,
        args: Self::Args,
    ) -> Result<Self, Self::Error>;
}

pub trait ArgValidate {
    fn validate(&self);
}
