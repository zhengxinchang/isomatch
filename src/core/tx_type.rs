use core::fmt;

#[derive(Copy, Clone, Debug)]
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
