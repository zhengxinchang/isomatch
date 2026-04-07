use crate::core::tx_base_error::TxBaseError;

/// Flags for TxBase.
/// bit 0: strand (0 for +, 1 for -)
/// bit 1: seq_hash is valid
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
#[repr(transparent)]
pub struct TxBaseFlags(pub u16);

impl TxBaseFlags {
    const NEG_STRAND_BIT: u16 = 1 << 0;
    const HAS_SEQ_HASH_BIT: u16 = 1 << 1;

    pub fn new(strand: u8, seq_has_hash: bool) -> Result<Self, TxBaseError> {
        let mut flags = Self(0);
        match strand {
            0 => {}
            1 => flags.0 |= Self::NEG_STRAND_BIT,
            _ => return Err(TxBaseError::InvalidStrand { strand }),
        }

        if seq_has_hash {
            flags.0 |= Self::HAS_SEQ_HASH_BIT;
        }

        Ok(flags)
    }

    pub fn get_strand(self) -> u8 {
        if self.0 & Self::NEG_STRAND_BIT == Self::NEG_STRAND_BIT {
            1
        } else {
            0
        }
    }

    pub fn set_strand(&mut self, strand: u8) -> Result<(), TxBaseError> {
        match strand {
            0 => {
                self.0 &= !Self::NEG_STRAND_BIT;
                Ok(())
            }
            1 => {
                self.0 |= Self::NEG_STRAND_BIT;
                Ok(())
            }
            _ => Err(TxBaseError::InvalidStrand { strand }),
        }
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

    pub fn strand(self) -> u8 {
        self.get_strand()
    }

    pub fn seq_has_hash(self) -> bool {
        self.get_seq_has_hash()
    }
}
