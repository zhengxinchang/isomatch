use flate2::bufread::MultiGzDecoder;
use serde::Serialize;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use xxhash_rust::xxh3::xxh3_128;
pub fn greetings2<T: Serialize>(msg: &T) {
    match serde_json::to_string_pretty(&msg) {
        Ok(json) => eprintln!("Parsed arguments:\n{}", json),
        Err(e) => eprintln!("Failed to print arguments: {}", e),
    }
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
    let mut file_reader = BufReader::new(File::open(path)?);
    let is_gzip = file_reader.fill_buf()?.starts_with(&[0x1f, 0x8b]);

    if is_gzip {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file_reader))))
    } else {
        Ok(Box::new(file_reader))
    }
}
