mod bucket;
mod reader;
use crate::error::Error;
use crate::error::Result;
use core::fmt;

pub use reader::{ChromID, GtfProfile, GtfReader, Transcript};

#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Strand {
    Plus = 0,
    Minus = 1,
    Unknown = 2,
}

impl Strand {
    pub fn to_bit(self) -> u8 {
        self as u8
    }

    pub fn from_bit(value: u8) -> Result<Self> {
        match value {
            0 => Ok(Strand::Plus),
            1 => Ok(Strand::Minus),
            2 => Ok(Strand::Unknown),
            _ => Err(Error::InvalidStrand { strand: value }),
        }
    }
}

impl From<Strand> for u8 {
    fn from(value: Strand) -> Self {
        value as u8
    }
}

impl TryFrom<u8> for Strand {
    type Error = Error;

    fn try_from(value: u8) -> Result<Self> {
        Self::from_bit(value)
    }
}

impl From<Strand> for char {
    fn from(s: Strand) -> Self {
        match s {
            Strand::Plus => '+',
            Strand::Minus => '-',
            Strand::Unknown => '.',
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Plus => write!(f, "Plus Strand"),
            Strand::Minus => write!(f, "Minus Strand"),
            Strand::Unknown => write!(f, "Unknown Strand"),
        }
    }
}

// valid record line from GTF input
pub struct GtfRecord<'a> {
    pub seqid: &'a str,
    pub feature: &'a str,
    pub start: u32,
    pub end: u32,
    pub strand: Strand,
    pub attributes: &'a str,
}

pub fn parse_record(line: &str) -> Result<GtfRecord<'_>> {
    let mut cols = line.splitn(9, '\t');
    let mut next = |field| {
        cols.next().ok_or_else(|| Error::InvalidGtfField {
            field,
            reason: "missing column".to_string(),
        })
    };

    let seqid = next(0)?;
    let _source = next(1)?;
    let feature = next(2)?;
    let start = next(3)?
        .parse::<u32>()
        .map_err(|_| Error::InvalidGtfField {
            field: 3,
            reason: "invalid start coordinate".to_string(),
        })?;
    let end = next(4)?
        .parse::<u32>()
        .map_err(|_| Error::InvalidGtfField {
            field: 4,
            reason: "invalid end coordinate".to_string(),
        })?;
    let _score = next(5)?;
    let strand = match next(6)? {
        "+" => Strand::Plus,
        "-" => Strand::Minus,
        "." => Strand::Unknown,
        _ => {
            return Err(Error::InvalidGtfField {
                field: 6,
                reason: "invalid strand".to_string(),
            });
        }
    };
    let _frame = next(7)?;

    let attributes = next(8)?.split('\t').next().unwrap_or("");

    // fix it after confirm isom
    // let attributes = next(8)?;

    Ok(GtfRecord {
        seqid,
        feature,
        start,
        end,
        strand,
        attributes,
    })
}

pub fn attribute<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    let mut saw_empty_match = false;

    for attr in attrs.split(';') {
        let attr = attr.trim();

        let mut parts = attr.splitn(2, char::is_whitespace);

        if parts.next() != Some(key) {
            continue;
        }

        let quoted = attr.find('"').and_then(|start| {
            attr[start + 1..]
                .find('"')
                .map(|len| &attr[start + 1..start + 1 + len])
        });
        let value = quoted.unwrap_or_else(|| attr.split_ascii_whitespace().nth(1).unwrap_or(""));

        if !value.is_empty() {
            return Some(value);
        }
        saw_empty_match = true;
    }

    saw_empty_match.then_some("")
}

#[cfg(test)]
mod tests {
    use super::{Strand, attribute, parse_record};
    use crate::error::Error;

    #[test]
    fn parse_record_reads_gtf_fields() {
        let line = "chr1\tsrc\texon\t12\t34\t.\t-\t.\tgene_id \"G1\";";
        let record = parse_record(line).unwrap();
        assert_eq!(record.seqid, "chr1");
        assert_eq!(record.feature, "exon");
        assert_eq!((record.start, record.end), (12, 34));
        assert_eq!(record.strand, Strand::Minus);
        assert_eq!(record.attributes, "gene_id \"G1\";");
    }

    #[test]
    fn parse_record_reports_invalid_fields() {
        for (line, expected_field) in [
            ("chr1\tsrc\texon", 3),
            ("chr1\tsrc\texon\tx\t34\t.\t+\t.\tgene_id G1;", 3),
            ("chr1\tsrc\texon\t12\tx\t.\t+\t.\tgene_id G1;", 4),
            ("chr1\tsrc\texon\t12\t34\t.\t?\t.\tgene_id G1;", 6),
        ] {
            match parse_record(line) {
                Err(Error::InvalidGtfField { field, .. }) => assert_eq!(field, expected_field),
                _ => panic!("expected invalid GTF record: {line}"),
            }
        }
    }

    #[test]
    fn attribute_preserves_legacy_lookup_behavior() {
        let cases = [
            ("gene_id_version \"V\"; gene_id \"G1\";", Some("G1")),
            ("gene_id G1;", Some("G1")),
            ("gene_id \"A B\";", Some("A B")),
            ("gene_id \"\"; gene_id \"G1\";", Some("G1")),
            ("gene_id \"\"; gene_id;", Some("")),
            ("gene_id_version \"V\";", None),
        ];

        for (attrs, expected) in cases {
            assert_eq!(attribute(attrs, "gene_id"), expected, "{attrs}");
        }
    }
}
