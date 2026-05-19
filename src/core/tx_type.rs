use core::fmt;

#[derive(Copy, Clone, Debug)]
pub enum TxType {
    MONO,
    ALLC,
    PRTC,
    NOTC,
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
