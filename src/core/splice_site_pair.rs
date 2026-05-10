use crate::{
    core::{core_error::TxBaseError, tx_strand::ISOMSTRAND},
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
            ISOMSTRAND::Unknown => {
                let plus_pack = Self::pack_for_known_strand(left, right, ISOMSTRAND::Plus);
                if Self::packed_is_canonical(plus_pack) {
                    return Ok(Self(plus_pack));
                }

                let minus_pack = Self::pack_for_known_strand(left, right, ISOMSTRAND::Minus);
                if Self::packed_is_canonical(minus_pack) {
                    return Ok(Self(minus_pack));
                }

                Ok(Self(5))
            }
            _ => Ok(Self(Self::pack_for_known_strand(left, right, strand))),
        }
    }

    pub fn from_packed(p: u8) -> Result<Self, TxBaseError> {
        Ok(Self(p))
    }

    pub fn is_canonical(&self) -> bool {
        Self::packed_is_canonical(self.0)
    }

    fn pack_for_known_strand(left: &[u8], right: &[u8], strand: ISOMSTRAND) -> u8 {
        let norm_left = normalized_site(left, &strand);
        let norm_right = normalized_site(right, &strand);

        let left_code = Self::site_code(&norm_left[0..2]);
        let right_code = Self::site_code(&norm_right[0..2]);

        if strand == ISOMSTRAND::Minus {
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
        (donor == 0 && acceptor == 1)
            || (donor == 2 && acceptor == 1)
            || (donor == 3 && acceptor == 4)
    }
}
