use xxhash_rust::xxh3::xxh3_128;

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

pub fn upper_nuc(b: u8) -> u8 {
    match b {
        b'a' => b'A',
        b't' => b'T',
        b'c' => b'C',
        b'g' => b'G',
        b'n' => b'N',
        other => other,
    }
}

pub fn hash_u8_slice(v: &[u8]) -> u128 {
    xxh3_128(v)
}
