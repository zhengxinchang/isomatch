use std::{
    fs::File,
    io::{self, BufRead, BufReader, ErrorKind, Read, Write},
    path::Path,
};

use crate::gtf::Strand;
use flate2::read::MultiGzDecoder;

const BUFREADER_CAPACITY: usize = 128 * 1024;

pub fn write_u8<W: Write>(writer: &mut W, value: u8) -> io::Result<()> {
    writer.write_all(&[value])
}

pub fn write_u16<W: Write>(writer: &mut W, value: u16) -> io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

pub fn write_u32<W: Write>(writer: &mut W, value: u32) -> io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

pub fn write_u64<W: Write>(writer: &mut W, value: u64) -> io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

pub fn write_bytes<W: Write>(writer: &mut W, bytes: &[u8]) -> io::Result<()> {
    let len = u32::try_from(bytes.len()).map_err(|_| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("field length {} exceeded u32", bytes.len()),
        )
    })?;
    write_u32(writer, len)?;
    writer.write_all(bytes)
}

pub fn read_u8<R: Read>(reader: &mut R) -> io::Result<u8> {
    let mut buf = [0u8; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

pub fn read_u8_opt<R: Read>(reader: &mut R) -> io::Result<Option<u8>> {
    let mut buf = [0u8; 1];
    match reader.read_exact(&mut buf) {
        Ok(()) => Ok(Some(buf[0])),
        Err(err) if err.kind() == ErrorKind::UnexpectedEof => Ok(None),
        Err(err) => Err(err),
    }
}

pub fn read_u16<R: Read>(reader: &mut R) -> io::Result<u16> {
    let mut buf = [0u8; 2];
    reader.read_exact(&mut buf)?;
    Ok(u16::from_le_bytes(buf))
}

pub fn read_u32<R: Read>(reader: &mut R) -> io::Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

pub fn read_u32_opt<R: Read>(reader: &mut R) -> io::Result<Option<u32>> {
    let mut buf = [0u8; 4];
    match reader.read_exact(&mut buf) {
        Ok(()) => Ok(Some(u32::from_le_bytes(buf))),
        Err(err) if err.kind() == ErrorKind::UnexpectedEof => Ok(None),
        Err(err) => Err(err),
    }
}

pub fn read_u64<R: Read>(reader: &mut R) -> io::Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

pub fn read_bytes<R: Read>(reader: &mut R) -> io::Result<Vec<u8>> {
    let len = read_u32(reader)? as usize;
    let mut bytes = vec![0u8; len];
    reader.read_exact(&mut bytes)?;
    Ok(bytes)
}

pub fn read_strand<R: Read>(reader: &mut R) -> io::Result<Strand> {
    Strand::try_from(read_u8(reader)?)
        .map_err(|err| io::Error::new(ErrorKind::InvalidData, err.to_string()))
}

pub fn bytes_to_string(bytes: Vec<u8>) -> io::Result<String> {
    String::from_utf8(bytes).map_err(|err| io::Error::new(ErrorKind::InvalidData, err.to_string()))
}

pub fn bytes_to_string_io(bytes: Vec<u8>) -> io::Result<String> {
    String::from_utf8(bytes).map_err(|err| io::Error::new(ErrorKind::InvalidData, err.to_string()))
}

pub fn open_file_bufread<P: AsRef<Path>>(path: P) -> std::io::Result<Box<dyn BufRead>> {
    let mut file_reader = BufReader::with_capacity(BUFREADER_CAPACITY, File::open(path)?);
    let is_gzip = file_reader.fill_buf()?.starts_with(&[0x1f, 0x8b]);

    if is_gzip {
        Ok(Box::new(BufReader::with_capacity(
            BUFREADER_CAPACITY,
            MultiGzDecoder::new(file_reader),
        )))
    } else {
        Ok(Box::new(file_reader))
    }
}
