use core::fmt;

#[derive(Copy, Clone, Debug)]
#[repr(u8)]
pub enum TxType {
    MONO,
    ALLC,
    PRTC,
    NOTC,
}

impl TxType {
    pub fn from_str(s: &str) -> Option<TxType> {
        match s {
            "MONO" => Some(TxType::MONO),
            "ALL_CA" => Some(TxType::ALLC),
            "PRT_CA" => Some(TxType::PRTC),
            "NOT_CA" => Some(TxType::NOTC),
            _ => None,
        }
    }

    pub fn to_u8(self) -> u8 {
        self as u8
    }

    pub fn from_u8(v: u8) -> Option<TxType> {
        match v {
            0 => Some(TxType::MONO),
            1 => Some(TxType::ALLC),
            2 => Some(TxType::PRTC),
            3 => Some(TxType::NOTC),
            _ => None,
        }
    }
}

impl fmt::Display for TxType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            TxType::MONO => "MONO",
            TxType::ALLC => "ALL_CA",
            TxType::PRTC => "PRT_CA",
            TxType::NOTC => "NOT_CA",
        };
        write!(f, "{s}")
    }
}
