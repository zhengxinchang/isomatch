use std::io::{Read, Seek, Write};

pub trait DiskSize {
    const DISK_SIZE: usize;
}
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

pub trait Decodable: Sized {
    type Error;
    type Args;
    fn decode_from<R: Read + Seek>(reader: &mut R, args: Self::Args) -> Result<Self, Self::Error>;
}

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

pub trait LogMemSize {
    fn get_mem_size(&self) -> usize;
}
