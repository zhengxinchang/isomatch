use std::{
    fs::{self, File},
    io::{self, BufReader, BufWriter, Cursor, ErrorKind, Read, Write},
    path::{Path, PathBuf},
};

use crate::gtf::Strand;
use crate::gtf::reader::Transcript;

pub type ChromID = u16;

const REC_TX: u8 = 1;
const REC_EXON: u8 = 2;

use crate::io::*;

pub struct Bucket {
    writer: Option<BufWriter<File>>,
    reader: Option<BufReader<File>>,
}

impl Bucket {
    pub fn init_writer<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self {
            writer: Some(BufWriter::new(File::create(path.as_ref())?)),
            reader: None,
        })
    }

    pub fn init_reader<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self {
            writer: None,
            reader: Some(BufReader::new(File::open(path.as_ref())?)),
        })
    }

    pub fn flush(&mut self) -> io::Result<()> {
        if let Some(writer) = self.writer.as_mut() {
            writer.flush()?;
        }
        Ok(())
    }

    pub fn write_tx(&mut self, tx: &TmpTxRec) -> io::Result<()> {
        let mut payload = Vec::new();
        write_u64(&mut payload, tx.hash)?;
        write_bytes(&mut payload, &tx.chrom)?;
        write_u32(&mut payload, tx.start)?;
        write_u32(&mut payload, tx.end)?;
        write_u8(&mut payload, tx.strand as u8)?;
        write_bytes(&mut payload, &tx.tx_id)?;
        write_bytes(&mut payload, &tx.gene_id)?;
        write_bytes(&mut payload, &tx.raw_attr_string)?;
        self.write_record(REC_TX, &payload)
    }

    pub fn write_exon(&mut self, exon: &TmpExonRec) -> io::Result<()> {
        let mut payload = Vec::new();
        write_u64(&mut payload, exon.hash)?;
        write_bytes(&mut payload, &exon.chrom)?;
        write_u32(&mut payload, exon.start)?;
        write_u32(&mut payload, exon.end)?;
        write_u8(&mut payload, exon.strand as u8)?;
        write_bytes(&mut payload, &exon.tx_id)?;
        write_bytes(&mut payload, &exon.gene_id)?;
        self.write_record(REC_EXON, &payload)
    }

    pub fn read_one(&mut self) -> io::Result<Option<TmpRec>> {
        let Some(reader) = self.reader.as_mut() else {
            return Err(io::Error::new(
                ErrorKind::InvalidInput,
                "bucket is not open for reading",
            ));
        };

        let Some(kind) = read_u8_opt(reader)? else {
            return Ok(None);
        };
        let payload_len = read_u32(reader)? as usize;
        let mut payload = vec![0u8; payload_len];
        reader.read_exact(&mut payload)?;
        let mut cursor = Cursor::new(payload);

        match kind {
            REC_TX => Ok(Some(TmpRec::Tx(TmpTxRec {
                hash: read_u64(&mut cursor)?,
                chrom: read_bytes(&mut cursor)?,
                start: read_u32(&mut cursor)?,
                end: read_u32(&mut cursor)?,
                strand: read_strand(&mut cursor)?,
                tx_id: read_bytes(&mut cursor)?,
                gene_id: read_bytes(&mut cursor)?,
                raw_attr_string: read_bytes(&mut cursor)?,
            }))),
            REC_EXON => Ok(Some(TmpRec::Exon(TmpExonRec {
                hash: read_u64(&mut cursor)?,
                chrom: read_bytes(&mut cursor)?,
                start: read_u32(&mut cursor)?,
                end: read_u32(&mut cursor)?,
                strand: read_strand(&mut cursor)?,
                tx_id: read_bytes(&mut cursor)?,
                gene_id: read_bytes(&mut cursor)?,
            }))),
            _ => Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("unknown bucket record kind {kind}"),
            )),
        }
    }

    fn write_record(&mut self, kind: u8, payload: &[u8]) -> io::Result<()> {
        let Some(writer) = self.writer.as_mut() else {
            return Err(io::Error::new(
                ErrorKind::InvalidInput,
                "bucket is not open for writing",
            ));
        };
        write_u8(writer, kind)?;
        write_bytes(writer, payload)
    }
}

pub enum TmpRec {
    Tx(TmpTxRec),
    Exon(TmpExonRec),
}

pub struct TmpTxRec {
    pub hash: u64,
    pub chrom: Vec<u8>,
    pub start: u32,
    pub end: u32,
    pub strand: Strand,
    pub tx_id: Vec<u8>,
    pub gene_id: Vec<u8>,
    pub raw_attr_string: Vec<u8>,
}

pub struct TmpExonRec {
    pub hash: u64,
    pub chrom: Vec<u8>,
    pub start: u32,
    pub end: u32,
    pub strand: Strand,
    pub tx_id: Vec<u8>,
    pub gene_id: Vec<u8>,
}

pub struct SortedBucket {
    pub writer: Option<BufWriter<File>>,
    pub reader: Option<BufReader<File>>,
}

impl SortedBucket {
    pub fn init_writer<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self {
            writer: Some(BufWriter::new(File::create(path.as_ref())?)),
            reader: None,
        })
    }

    pub fn init_reader<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self {
            writer: None,
            reader: Some(BufReader::new(File::open(path.as_ref())?)),
        })
    }

    pub fn dump_tx_structure(&mut self, tx: &Transcript) -> io::Result<()> {
        let Some(writer) = self.writer.as_mut() else {
            return Err(io::Error::new(
                ErrorKind::InvalidInput,
                "sorted bucket is not open for writing",
            ));
        };

        let mut payload = Vec::new();
        write_u16(&mut payload, tx.chrom_id)?;
        write_u32(&mut payload, tx.start)?;
        write_u32(&mut payload, tx.end)?;
        write_u8(&mut payload, tx.strand as u8)?;
        write_bytes(&mut payload, tx.tx_id.as_bytes())?;
        write_bytes(&mut payload, tx.gene_id.as_bytes())?;
        write_u32(&mut payload, tx.exons.len() as u32)?;
        for (start, end) in &tx.exons {
            write_u32(&mut payload, *start)?;
            write_u32(&mut payload, *end)?;
        }
        match &tx.attr_string {
            Some(attr) => {
                write_u8(&mut payload, 1)?;
                write_bytes(&mut payload, attr)?;
            }
            None => write_u8(&mut payload, 0)?,
        }
        write_bytes(writer, &payload)
    }

    pub fn flush(&mut self) -> io::Result<()> {
        if let Some(writer) = self.writer.as_mut() {
            writer.flush()?;
        }
        Ok(())
    }

    pub fn read_one(&mut self) -> io::Result<Option<Transcript>> {
        let Some(reader) = self.reader.as_mut() else {
            return Err(io::Error::new(
                ErrorKind::InvalidInput,
                "sorted bucket is not open for reading",
            ));
        };

        let Some(payload_len) = read_u32_opt(reader)? else {
            return Ok(None);
        };
        let mut payload = vec![0u8; payload_len as usize];
        reader.read_exact(&mut payload)?;
        let mut cursor = Cursor::new(payload);

        let chrom_id = read_u16(&mut cursor)?;
        let start = read_u32(&mut cursor)?;
        let end = read_u32(&mut cursor)?;
        let strand = read_strand(&mut cursor)?;
        let tx_id = bytes_to_string_io(read_bytes(&mut cursor)?)?;
        let gene_id = bytes_to_string_io(read_bytes(&mut cursor)?)?;
        let exon_count = read_u32(&mut cursor)? as usize;
        let mut exons = Vec::with_capacity(exon_count);
        for _ in 0..exon_count {
            exons.push((read_u32(&mut cursor)?, read_u32(&mut cursor)?));
        }
        let attr_string = match read_u8(&mut cursor)? {
            0 => None,
            1 => Some(read_bytes(&mut cursor)?),
            value => {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    format!("invalid attr presence flag {value}"),
                ));
            }
        };

        Ok(Some(Transcript {
            gidx: 0,
            chrom_id,
            start,
            end,
            strand,
            exons,
            tx_id,
            gene_id,
            is_empty: false,
            attr_string,
        }))
    }
}

/// intermediate structure for aggregating exon and transcript records in one TxStrcture
#[derive(Default)]
pub struct Rec2TxStrctureTmp {
    pub tx: Option<Transcript>,
    pub transcript_chrom_id: Option<ChromID>,
    pub attr_string: Option<Vec<u8>>,
}

pub struct TempDirCleanup {
    path: PathBuf,
    armed: bool,
}

impl TempDirCleanup {
    pub fn new(path: PathBuf) -> Self {
        Self { path, armed: true }
    }

    pub fn disarm(&mut self) {
        self.armed = false;
    }
}

impl Drop for TempDirCleanup {
    fn drop(&mut self) {
        if self.armed {
            let _ = fs::remove_dir_all(&self.path);
        }
    }
}
