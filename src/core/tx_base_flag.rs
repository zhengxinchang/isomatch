use crate::core::tx_strand::ISOMSTRAND;
use crate::core::tx_base_error::TxBaseError;

/// Flags for TxBase.
/// bit 0-1: strand (0 for +, 1 for -, 2 for unknown)
/// bit 2: seq_hash is valid
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
#[repr(transparent)]
pub struct TxBaseFlags(pub u16);

impl TxBaseFlags {
    const STRAND_MASK: u16 = 0b11;
    const HAS_SEQ_HASH_BIT: u16 = 1 << 2;

    pub fn new(strand: ISOMSTRAND, seq_has_hash: bool) -> Result<Self, TxBaseError> {
        let mut flags = Self(u16::from(strand.to_bit()));

        if seq_has_hash {
            flags.0 |= Self::HAS_SEQ_HASH_BIT;
        }

        Ok(flags)
    }

    pub fn get_strand(self) -> ISOMSTRAND {
        ISOMSTRAND::from_bit((self.0 & Self::STRAND_MASK) as u8).unwrap()
    }

    pub fn set_strand(&mut self, strand: ISOMSTRAND) -> Result<(), TxBaseError> {
        self.0 &= !Self::STRAND_MASK;
        self.0 |= u16::from(strand.to_bit());
        Ok(())
    }

    pub fn get_seq_has_hash(&self) -> bool {
        self.0 & Self::HAS_SEQ_HASH_BIT == Self::HAS_SEQ_HASH_BIT
    }

    pub fn set_seq_has_hash(&mut self, has_hash: bool) {
        if has_hash {
            self.0 |= Self::HAS_SEQ_HASH_BIT;
        } else {
            self.0 &= !Self::HAS_SEQ_HASH_BIT;
        }
    }

    pub fn bits(self) -> u16 {
        self.0
    }

    // pub fn strand(self) -> u8 {
    //     self.get_strand()
    // }

    pub fn seq_has_hash(self) -> bool {
        self.get_seq_has_hash()
    }
}
