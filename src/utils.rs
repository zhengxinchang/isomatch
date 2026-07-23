use ahash::RandomState;
use flate2::bufread::MultiGzDecoder;
use serde::Serialize;
use std::collections::HashMap;
use std::format;
use std::fs::File;
use std::hash::Hash;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::LazyLock;
use thiserror::Error;
use xxhash_rust::xxh3::xxh3_128;

use crate::constants;
use crate::{
    core::tx_strand::ISOMSTRAND,
    index::{attributes_index::AttrIndexReader, reader::IndexReader},
};

const BUFREADER_CAPACITY: usize = 128 * 1024;

#[derive(Error, Debug)]
pub enum UtilsError {
    #[error(transparent)]
    Io(#[from] std::io::Error),

    #[error("Missing {key} in ISOM sample header: {line}")]
    MissingSampleHeaderValue { key: String, line: String },

    #[error("Invalid ISOM sample id {id}: {reason}")]
    InvalidSampleId { id: String, reason: String },
}

#[derive(Clone, Debug)]
pub struct IsomSample {
    pub file_name: String,
}

pub fn print_json_block<T: Serialize>(title: &str, msg: &T) {
    match serde_json::to_string_pretty(&msg) {
        Ok(json) => eprintln!("{}:\n{}", title, json),
        Err(e) => eprintln!("Failed to print {}: {}", title, e),
    }
}

pub fn save_json_block<T: Serialize, P: AsRef<Path>>(path: &P, msg: &T) -> Result<(), io::Error> {
    let mut file = File::create(path)?;
    match serde_json::to_string_pretty(&msg) {
        Ok(json) => file.write_all(json.as_bytes()),
        Err(e) => Err(io::Error::new(
            io::ErrorKind::Other,
            format!("Failed to save {:?}: {}", path.as_ref(), e),
        )),
    }
}

pub fn greetings2<T: Serialize>(msg: &T) {
    print_json_block("Parsed arguments", msg);
}

pub fn require_file(
    label: &str,
    path: &Path,
    error_msg: &mut String,
    has_error: &mut bool,
) -> bool {
    if !path.exists() {
        error_msg.push_str(&format!("\n{} does not exist: {:?}", label, path));
        *has_error = true;
        return false;
    }

    if !path.is_file() {
        error_msg.push_str(&format!("\n{} is not a file: {:?}", label, path));
        *has_error = true;
        return false;
    }

    true
}

pub fn checksum_file(path: &Path) -> std::io::Result<([u8; 16], u64)> {
    let metadata = std::fs::metadata(path)?;
    let size = metadata.len();

    let mut file = std::fs::File::open(path)?;
    let mut buf = vec![0u8; 64 * 1024];
    let mut hasher = xxhash_rust::xxh3::Xxh3::new();
    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    let hash = hasher.digest128();

    Ok((hash.to_le_bytes(), size))
}

pub fn hash_str(s: &str) -> u128 {
    xxh3_128(s.as_bytes())
}

pub fn hash_u8_vec(v: &Vec<u8>) -> u128 {
    xxh3_128(v)
}

pub fn hash_u8_slice(v: &[u8]) -> u128 {
    xxh3_128(v)
}

static HASHER: LazyLock<RandomState> =
    LazyLock::new(|| RandomState::with_seeds(9336, 5920, 6784, 4496));

/// ONLY used in memory, DO NOT used in persistance purpose.
/// the hash is not stable in different paltform, language.
pub fn ahash_vec<T: Hash>(value: &T) -> u64 {
    HASHER.hash_one(value)
}

pub fn trim_chr_prefix_to_upper(chrom: &str) -> String {
    chrom
        .to_ascii_uppercase()
        .trim_start_matches("CHR")
        .to_string()
}

pub fn pad_chrom_prefix(chrom: &str) -> String {
    if chrom.starts_with("chr") {
        chrom.to_string()
    } else {
        format!("chr{}", chrom)
    }
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

pub fn rev_comp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

pub fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' => b'N',
        other => other,
    }
}

fn upper_nuc(b: u8) -> u8 {
    match b {
        b'a' => b'A',
        b't' => b'T',
        b'c' => b'C',
        b'g' => b'G',
        b'n' => b'N',
        other => other,
    }
}
/// reverse site acoording to strand
/// also convert bases to upaer cases
pub fn normalized_site(site: &[u8], strand: &ISOMSTRAND) -> Vec<u8> {
    match strand {
        ISOMSTRAND::Minus => rev_comp(site),
        _ => site.iter().map(|&b| upper_nuc(b)).collect(),
    }
}

pub fn is_gzipped(p: &Path) -> bool {
    p.extension()
        .and_then(|s| s.to_str())
        .map(|s| s.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
}

pub fn is_isomatch_merged_gtf<P: AsRef<Path>>(path: P) -> Result<bool, UtilsError> {
    let mut reader = open_file_bufread(path)?;
    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        if !line.starts_with('#') {
            break;
        }
        if line.starts_with("##ISOM <VERSION>")
            && line.contains(&format!("schema=\"{}\"", constants::ISOM_GTF_SCHEMA))
        {
            return Ok(true);
        }
        line.clear();
    }
    Ok(false)
}

pub fn read_samples<P: AsRef<Path>>(path: P) -> Result<HashMap<u32, IsomSample>, UtilsError> {
    let mut samples = HashMap::new();
    let mut reader = open_file_bufread(path)?;
    let mut line = String::new();

    while reader.read_line(&mut line)? != 0 {
        if !line.starts_with('#') {
            break;
        }

        if line.starts_with("##ISOM <SAMPLE>") {
            let id = parse_header_value(&line, "id").ok_or_else(|| {
                UtilsError::MissingSampleHeaderValue {
                    key: "id".to_string(),
                    line: line.trim_end().to_string(),
                }
            })?;
            let input = parse_header_value(&line, "input").ok_or_else(|| {
                UtilsError::MissingSampleHeaderValue {
                    key: "input".to_string(),
                    line: line.trim_end().to_string(),
                }
            })?;
            let file_id = parse_sample_id(&id)?;
            let file_name = PathBuf::from(&input)
                .file_name()
                .map(|name| name.to_string_lossy().to_string())
                .unwrap_or(input);
            samples.insert(file_id, IsomSample { file_name });
        }

        line.clear();
    }

    Ok(samples)
}

fn parse_sample_id(id: &str) -> Result<u32, UtilsError> {
    match id.trim_start_matches('S').parse::<u32>() {
        Ok(0) => Err(UtilsError::InvalidSampleId {
            id: id.to_string(),
            reason: "sample id 0 is invalid".to_string(),
        }),
        Ok(file_id) => Ok(file_id),
        Err(err) => Err(UtilsError::InvalidSampleId {
            id: id.to_string(),
            reason: err.to_string(),
        }),
    }
}

fn parse_header_value(line: &str, key: &str) -> Option<String> {
    let needle = format!("{key}=\"");
    let start = line.find(&needle)? + needle.len();
    let end = line[start..].find('"')? + start;
    Some(line[start..end].replace("\\\"", "\"").replace("\\\\", "\\"))
}

pub fn check_index_ready<P: AsRef<Path>>(gtf_path: P) -> bool {
    let gtf_path = gtf_path.as_ref();

    let Ok(metadata) = std::fs::metadata(gtf_path) else {
        return false;
    };

    let mut isomx_path = gtf_path.to_path_buf();
    isomx_path.add_extension("isomx");
    let index_header = match IndexReader::load_header(isomx_path) {
        Ok(header) => header,
        Err(_) => return false,
    };
    if index_header.version != constants::ISOMX_VERSION
        || index_header.gtf_file_size != metadata.len()
    {
        return false;
    }

    let mut isoms_path = gtf_path.to_path_buf();
    isoms_path.add_extension("isoms");
    let attr_header = match AttrIndexReader::load_header(isoms_path) {
        Ok(header) => header,
        Err(_) => return false,
    };
    if attr_header.version != constants::ISOMS_VERSION || attr_header.md5 != index_header.md5 {
        return false;
    }

    let mut reader = match open_file_bufread(gtf_path) {
        Ok(reader) => reader,
        Err(_) => return false,
    };
    let mut hasher = xxhash_rust::xxh3::Xxh3::new();
    let mut buf = vec![0u8; 64 * 1024];
    loop {
        let n = match reader.read(&mut buf) {
            Ok(n) => n,
            Err(_) => return false,
        };
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }

    index_header.md5 == hasher.digest128().to_le_bytes()
}
