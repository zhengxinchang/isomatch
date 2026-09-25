use crate::{
    traits::{Encodable, PartialLoad},
    utils::rev_comp,
    utils::upper_nuc,
};

use crate::gtf::Strand;

use crate::error::IndexDataError;

#[derive(Debug, Clone, Default)]
pub struct SpliceSitePool {
    sites: Vec<SpliceSitePair>,
}

impl SpliceSitePool {
    pub fn new() -> Self {
        Self { sites: Vec::new() }
    }

    pub fn with_capacity(capacity: usize) -> Result<Self, IndexDataError> {
        Ok(Self {
            sites: Vec::with_capacity(capacity),
        })
    }

    pub fn add_pairs(
        &mut self,
        pairs: &[(Vec<u8>, Vec<u8>)],
        strand: Strand,
    ) -> Result<SpliceSiteSpan, IndexDataError> {
        let offset = u32::try_from(self.sites.len()).map_err(|_| IndexDataError::PoolTooLarge)?;
        let count = u16::try_from(pairs.len()).map_err(|_| IndexDataError::InvalidEncoding {
            msg: format!("too many splice site pairs: {}", pairs.len()),
        })?;

        for (left_site, right_site) in pairs {
            self.sites
                .push(SpliceSitePair::pack(left_site, right_site, strand)?);
        }

        Ok(SpliceSiteSpan { offset, count })
    }

    pub fn get_pair(&self, span: SpliceSiteSpan) -> Result<&[SpliceSitePair], IndexDataError> {
        let start = usize::try_from(span.offset).map_err(|_| IndexDataError::InvalidSpan {
            offset: span.offset,
            count: span.count,
            pool_len: self.sites.len(),
        })?;
        let end = start + usize::from(span.count);

        self.sites
            .get(start..end)
            .ok_or(IndexDataError::InvalidSpan {
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
    type Error = IndexDataError;

    fn encode_to<W: std::io::Write>(&self, writer: &mut W) -> Result<usize, Self::Error> {
        let bytes: Vec<u8> = self.sites.iter().map(|pair| pair.0).collect();
        writer
            .write_all(&bytes)
            .map_err(|e| IndexDataError::Io(e.to_string()))?;
        Ok(bytes.len())
    }
}

impl PartialLoad for SpliceSitePool {
    type Error = IndexDataError;
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
            .map_err(IndexDataError::io)?;
        reader
            .read_exact(&mut buf)
            .map_err(|e| IndexDataError::Io(e.to_string()))?;

        let mut sites = Vec::with_capacity(buf.len());
        for byte in buf {
            sites.push(SpliceSitePair::from_packed(byte)?);
        }

        Ok(Self { sites })
    }
}

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
    pub fn pack(left: &[u8], right: &[u8], strand: Strand) -> Result<Self, IndexDataError> {
        if left.len() != 2 || right.len() != 2 {
            return Err(IndexDataError::InvalidSpliceSite {
                site: format!("{:?},{:?}", left, right),
            });
        }

        match strand {
            Strand::Unknown => {
                let plus_pack = Self::pack_for_known_strand(left, right, Strand::Plus);
                if Self::packed_is_canonical(plus_pack) {
                    return Ok(Self(plus_pack));
                }

                let minus_pack = Self::pack_for_known_strand(left, right, Strand::Minus);
                if Self::packed_is_canonical(minus_pack) {
                    return Ok(Self(minus_pack));
                }

                Ok(Self(5))
            }
            _ => Ok(Self(Self::pack_for_known_strand(left, right, strand))),
        }
    }

    pub fn from_packed(p: u8) -> Result<Self, IndexDataError> {
        Ok(Self(p))
    }

    pub fn is_canonical(&self) -> bool {
        Self::packed_is_canonical(self.0)
    }

    fn pack_for_known_strand(left: &[u8], right: &[u8], strand: Strand) -> u8 {
        let norm_left = normalized_site(left, &strand);
        let norm_right = normalized_site(right, &strand);

        let left_code = Self::site_code(&norm_left[0..2]);
        let right_code = Self::site_code(&norm_right[0..2]);

        if strand == Strand::Minus {
            right_code << 4 | left_code
        } else {
            left_code << 4 | right_code
        }
    }

    fn site_code(site: &[u8]) -> u8 {
        match site {
            [b'G', b'T'] => 0,
            [b'A', b'G'] => 1,
            [b'G', b'C'] => 2,
            [b'A', b'T'] => 3,
            [b'A', b'C'] => 4,
            _ => 5,
        }
    }

    fn packed_is_canonical(packed: u8) -> bool {
        let donor = packed >> 4;
        let acceptor = packed & 0x0F;
        (donor == 0 || donor == 2) && acceptor == 1 || (donor == 3 && acceptor == 4)
    }
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

/// reverse site acoording to strand
/// also convert bases to upaer cases
pub fn normalized_site(site: &[u8], strand: &Strand) -> Vec<u8> {
    match strand {
        Strand::Minus => rev_comp(site),
        _ => site.iter().map(|&b| upper_nuc(b)).collect(),
    }
}
