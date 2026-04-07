use std::fmt::Display;

pub enum TxType {
    MONO,
    ALLC,
    PRTC,
    NOTC,
}

pub enum MergeStatus {
    SNGL, // not merged, single
    REPR, // representive tx in the merging group
    EXAC, // exact match, including the junciton and TSS
    ENDS, // junction matches and TSS TES in wobble region
    M2CA, // non-canonical to canonical by junction wobble compare
    SXSF, // small exon boundary shift
    SXMS, // small exon mssing
}

impl Display for MergeStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let out = match self {
            MergeStatus::SNGL => &"SNGL",
            MergeStatus::REPR => &"REPR",
            MergeStatus::EXAC => &"EXAC",
            MergeStatus::ENDS => &"ENDS",
            MergeStatus::M2CA => &"M2CA",
            MergeStatus::SXSF => &"SXSF",
            MergeStatus::SXMS => &"SXMS",
        };
        write!(f, "{}", out)
    }
}
