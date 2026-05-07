use crate::{
    core::{tx_base_error::TxBaseError, tx_strand::ISOMSTRAND},
    traits::{Encodable, PartialLoad},
    utils::normalized_site,
};

/// Packed Splice Site
/// Negative strand bases will be reverse complement
/// Site projection:
/// GT --> 0
/// AG --> 1
/// GC --> 2
/// AT --> 3
/// AC --> 4
/// other -->5
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
#[repr(transparent)]
pub struct SpliceSitePair(pub u8);

impl SpliceSitePair {
    pub fn pack(left: &[u8], right: &[u8], strand: ISOMSTRAND) -> Result<Self, TxBaseError> {
        if left.len() != 2 || right.len() != 2 {
            return Err(TxBaseError::InvalidSpliceSite {
                site: format!("{:?},{:?}", left, right),
            });
        }

        match strand {
            ISOMSTRAND::Unknown => { // for unkown strand transcript, the splice site is not valid. always return 5
                Ok(Self(5))
            }
            _ => {
                let norm_left = normalized_site(left, &strand);
                let norm_right = normalized_site(right,& strand);

                let left_code: u8 = match norm_left[0..2] {
                    [b'G', b'T'] => 0,
                    [b'A', b'G'] => 1,
                    [b'G', b'C'] => 2,
                    [b'A', b'T'] => 3,
                    [b'A', b'C'] => 4,
                    _ => 5,
                };

                let right_code: u8 = match norm_right[0..2] {
                    [b'G', b'T'] => 0,
                    [b'A', b'G'] => 1,
                    [b'G', b'C'] => 2,
                    [b'A', b'T'] => 3,
                    [b'A', b'C'] => 4,
                    _ => 5,
                };

                let pack = if strand == ISOMSTRAND::Minus {
                    right_code << 4 | left_code
                } else {
                    left_code << 4 | right_code
                };
                Ok(Self(pack))
            }
        }
    }

    pub fn from_packed(p: u8) -> Result<Self, TxBaseError> {
        Ok(Self(p))
    }

    pub fn is_canonical(&self) -> bool {
        // canonical:
        // GT-AG --> 0 1
        // GC-AG --> 2 1
        // AT-AC --> 3 4
        // other 5
        let donor = self.0 >> 4;
        let acceptor = self.0 & 0x0F;
        (donor == 0 && acceptor == 1)
            || (donor == 2 && acceptor == 1)
            || (donor == 3 && acceptor == 4)
    }
    // /// This is
    // pub fn is_canonical_from_u8(c: u8) -> bool {
    //     let donor = c >> 4;
    //     let acceptor = c & 0x0F;
    //     (donor == 0 && acceptor == 1)
    //         || (donor == 2 && acceptor == 1)
    //         || (donor == 3 && acceptor == 4)
    // }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct SpliceSiteSpan {
    pub offset: u32,
    pub count: u16,
}

impl SpliceSiteSpan {
    pub fn is_empty(self) -> bool {
        self.count == 0
    }

    pub fn end_offset(self) -> u32 {
        self.offset + u32::from(self.count)
    }
}

#[derive(Debug, Clone, Default)]
pub struct SpliceSitePool {
    pub sites: Vec<SpliceSitePair>,
}

impl SpliceSitePool {
    pub fn new() -> Self {
        Self { sites: Vec::new() }
    }

    pub fn with_capacity(capacity: usize) -> Result<Self, TxBaseError> {
        Ok(Self {
            sites: Vec::with_capacity(capacity),
        })
    }

    pub fn add_pairs(
        &mut self,
        pairs: &[(Vec<u8>, Vec<u8>)],
        strand: ISOMSTRAND,
    ) -> Result<SpliceSiteSpan, TxBaseError> {
        let offset = u32::try_from(self.sites.len()).map_err(|_| TxBaseError::PoolTooLarge)?;
        let count = u16::try_from(pairs.len()).map_err(|_| TxBaseError::InvalidEncoding {
            msg: format!("too many splice site pairs: {}", pairs.len()),
        })?;

        for (left_site, right_site) in pairs {
            self.sites
                .push(SpliceSitePair::pack(&left_site, &right_site, strand)?);
        }

        Ok(SpliceSiteSpan { offset, count })
    }

    pub fn get_pair(&self, span: SpliceSiteSpan) -> Result<&[SpliceSitePair], TxBaseError> {
        let start = usize::try_from(span.offset).map_err(|_| TxBaseError::InvalidSpan {
            offset: span.offset,
            count: span.count,
            pool_len: self.sites.len(),
        })?;
        let end = start + usize::from(span.count);

        self.sites.get(start..end).ok_or(TxBaseError::InvalidSpan {
            offset: span.offset,
            count: span.count,
            pool_len: self.sites.len(),
        })
    }

    pub fn len(&self) -> usize {
        self.sites.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sites.is_empty()
    }
}

impl Encodable for SpliceSitePool {
    type Error = TxBaseError;

    fn encode_to<W: std::io::Write>(&self, writer: &mut W) -> Result<usize, Self::Error> {
        let bytes: Vec<u8> = self.sites.iter().map(|pair| pair.0).collect();
        writer
            .write_all(&bytes)
            .map_err(|e| TxBaseError::Io(e.to_string()))?;
        Ok(bytes.len())
    }
}

impl PartialLoad for SpliceSitePool {
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
            .map_err(TxBaseError::io)?;
        reader
            .read_exact(&mut buf)
            .map_err(|e| TxBaseError::Io(e.to_string()))?;

        let mut sites = Vec::with_capacity(buf.len());
        for byte in buf {
            sites.push(SpliceSitePair::from_packed(byte)?);
        }

        Ok(Self { sites })
    }
}

#[cfg(test)]
mod tests {
    use super::SpliceSitePair;
    use crate::core::tx_strand::ISOMSTRAND;

    #[test]
    fn canonical_pairs_on_plus_strand_are_recognized() {
        let gt_ag = SpliceSitePair::pack(b"GT", b"AG", ISOMSTRAND::Plus).unwrap();
        let gc_ag = SpliceSitePair::pack(b"GC", b"AG", ISOMSTRAND::Plus).unwrap();
        let at_ac = SpliceSitePair::pack(b"AT", b"AC", ISOMSTRAND::Plus).unwrap();

        assert!(gt_ag.is_canonical());
        assert!(gc_ag.is_canonical());
        assert!(at_ac.is_canonical());
    }

    #[test]
    fn canonical_pairs_on_minus_strand_are_recognized_after_normalization() {
        let gt_ag = SpliceSitePair::pack(b"CT", b"AC", ISOMSTRAND::Minus).unwrap();
        let gc_ag = SpliceSitePair::pack(b"CT", b"GC", ISOMSTRAND::Minus).unwrap();
        let at_ac = SpliceSitePair::pack(b"GT", b"AT", ISOMSTRAND::Minus).unwrap();

        assert!(gt_ag.is_canonical());
        assert!(gc_ag.is_canonical());
        assert!(at_ac.is_canonical());
    }

    #[test]
    fn noncanonical_pairs_are_not_marked_canonical() {
        let gt_gc = SpliceSitePair::pack(b"GT", b"GC", ISOMSTRAND::Plus).unwrap();
        let other = SpliceSitePair::pack(b"AA", b"AA", ISOMSTRAND::Plus).unwrap();

        assert!(!gt_gc.is_canonical());
        assert!(!other.is_canonical());
    }
}
