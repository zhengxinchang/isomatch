use core::fmt;

use crate::core::core_error::TxBaseError;

#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ISOMSTRAND {
    Plus = 0,
    Minus = 1,
    Unknown = 2,
}

impl ISOMSTRAND {
    pub fn to_bit(self) -> u8 {
        self as u8
    }

    pub fn from_bit(value: u8) -> Result<Self, TxBaseError> {
        match value {
            0 => Ok(ISOMSTRAND::Plus),
            1 => Ok(ISOMSTRAND::Minus),
            2 => Ok(ISOMSTRAND::Unknown),
            _ => Err(TxBaseError::InvalidStrand { strand: value }),
        }
    }
}

impl From<ISOMSTRAND> for u8 {
    fn from(value: ISOMSTRAND) -> Self {
        value as u8
    }
}

impl TryFrom<u8> for ISOMSTRAND {
    type Error = TxBaseError;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        Self::from_bit(value)
    }
}

impl From<ISOMSTRAND> for char {
    fn from(s: ISOMSTRAND) -> Self {
        match s {
            ISOMSTRAND::Plus => '+',
            ISOMSTRAND::Minus => '-',
            ISOMSTRAND::Unknown => '.',
        }
    }
}

impl fmt::Display for ISOMSTRAND {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ISOMSTRAND::Plus => write!(f, "Plus Strand"),
            ISOMSTRAND::Minus => write!(f, "Minus Strand"),
            ISOMSTRAND::Unknown => write!(f, "Unknown Strand"),
        }
    }
}
