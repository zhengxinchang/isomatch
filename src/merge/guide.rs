use crate::core::tx_strand::ISOMSTRAND;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

// bed format for guide tss tes
// chromosome      start   end     ID       score   strand
// chr1    16013   16020   rfhg_1.1        1       -

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GuideBEDType {
    Tss,
    Tes,
}

pub type ChromMap = FxHashMap<String, String>;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GuideInterval {
    pub start: u32, // 1-based closed
    pub end: u32,   // 1-based closed
    pub score: f32,
}

impl GuideInterval {
    #[inline]
    pub fn overlaps_point(&self, pos: u32) -> bool {
        self.start <= pos && pos <= self.end
    }

    #[inline]
    pub fn overlaps_range(&self, start: u32, end: u32) -> bool {
        self.start <= end && start <= self.end
    }

    #[inline]
    pub fn len(&self) -> u32 {
        self.end - self.start + 1
    }
}

#[derive(Debug, Clone, Default)]
pub struct ChromGuideIndex {
    intervals: Vec<GuideInterval>,
    max_len: u32,
}

impl ChromGuideIndex {
    pub fn intervals(&self) -> &[GuideInterval] {
        &self.intervals
    }

    pub fn max_len(&self) -> u32 {
        self.max_len
    }

    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }

    pub fn query_overlaps(&self, pos: u32) -> Vec<&GuideInterval> {
        self.query_overlaps_with_flank(pos, 0)
    }

    pub fn query_overlaps_with_flank(&self, pos: u32, flank: u32) -> Vec<&GuideInterval> {
        let start = pos.saturating_sub(flank);
        let end = pos.saturating_add(flank);
        self.query_overlaps_range(start, end)
    }

    fn query_overlaps_range(&self, start: u32, end: u32) -> Vec<&GuideInterval> {
        if self.intervals.is_empty() {
            return Vec::new();
        }

        let lower_start = start.saturating_sub(self.max_len);
        let lo = self
            .intervals
            .partition_point(|interval| interval.start < lower_start);
        let hi = self
            .intervals
            .partition_point(|interval| interval.start <= end);

        self.intervals[lo..hi]
            .iter()
            .filter(|interval| interval.overlaps_range(start, end))
            .collect()
    }
}

#[derive(Debug)]
pub struct GuideDb {
    guide_type: GuideBEDType,
    // bed_chroms: HashSet<String>,
    by_chrom_strand: FxHashMap<(String, ISOMSTRAND), ChromGuideIndex>,
    chrmap: Option<ChromMap>,
}

impl GuideDb {
    pub fn from_bed_path<P: AsRef<Path>>(
        path: P,
        guide_type: GuideBEDType,
        chrmap_path: &Option<P>,
    ) -> Result<Self, GuideError> {
        let path = path.as_ref();
        let file = File::open(path).map_err(|err| GuideError::Io {
            path: path.to_path_buf(),
            source: err,
        })?;
        let reader = BufReader::new(file);
        // Self::from_bed_reader(reader, guide_type)
        // let mut bed_chroms = HashSet::default();

        let mut grouped: FxHashMap<(String, ISOMSTRAND), Vec<GuideInterval>> = FxHashMap::default();

        for (line_no, line_result) in reader.lines().enumerate() {
            let raw_line = line_result.map_err(|err| GuideError::Io {
                path: path.to_path_buf(),
                source: err,
            })?;

            // Always skip the first line, which is expected to be the header.
            if line_no == 0 {
                continue;
            }

            let line = raw_line.trim();
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track")
                || line.starts_with("browser")
            {
                continue;
            }

            let record = parse_bed_record(line, line_no + 1)?;
            // bed_chroms.insert(record.chrom.clone());
            grouped
                .entry((record.chrom, record.strand))
                .or_default()
                .push(record.interval);
        }

        let by_chrom_strand = grouped
            .into_iter()
            .map(|(key, mut intervals)| {
                intervals.sort_by_key(|interval| interval.start);
                let max_len = intervals.iter().map(GuideInterval::len).max().unwrap_or(0);
                (key, ChromGuideIndex { intervals, max_len })
            })
            .collect();

        let chrmap = if let Some(p) = chrmap_path {
            Some(load_chrmap_path(p)?)
        } else {
            None
        };

        Ok(Self {
            guide_type,
            by_chrom_strand,
            // bed_chroms,
            chrmap,
        })
    }

    pub fn guide_type(&self) -> GuideBEDType {
        self.guide_type
    }

    pub fn get_index(&self, chrom: &str, strand: ISOMSTRAND) -> Option<&ChromGuideIndex> {
        if let Some(index) = self.by_chrom_strand.get(&(chrom.to_string(), strand)) {
            return Some(index);
        }

        if let Some(secondary_chrom) = self.chrmap.as_ref().and_then(|map| map.get(chrom)) {
            return self.by_chrom_strand.get(&(secondary_chrom.clone(), strand));
        }

        None
    }

    pub fn query_overlaps(&self, chrom: &str, strand: ISOMSTRAND, pos: u32) -> Vec<&GuideInterval> {
        self.get_index(chrom, strand)
            .map(|index| index.query_overlaps(pos))
            .unwrap_or_default()
    }

    pub fn query_overlaps_with_flank(
        &self,
        chrom: &str,
        strand: &ISOMSTRAND,
        pos: u32,
        flank: u32,
    ) -> Vec<&GuideInterval> {
        self.get_index(chrom, *strand)
            .map(|index| index.query_overlaps_with_flank(pos, flank))
            .unwrap_or_default()
    }

    pub fn len(&self) -> usize {
        self.by_chrom_strand.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_chrom_strand.is_empty()
    }

    pub fn chrmap(&self) -> Option<&ChromMap> {
        self.chrmap.as_ref()
    }
}

pub fn load_chrmap_path<P: AsRef<Path>>(path: P) -> Result<ChromMap, GuideError> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|err| GuideError::Io {
        path: path.to_path_buf(),
        source: err,
    })?;
    let reader = BufReader::new(file);
    let mut chrmap = FxHashMap::default();

    for (line_no, line_result) in reader.lines().enumerate() {
        let raw_line = line_result.map_err(|err| GuideError::Io {
            path: path.to_path_buf(),
            source: err,
        })?;

        // Always skip the first line, which is expected to be the header.
        if line_no == 0 {
            continue;
        }

        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let mut fields = line.split_whitespace();
        let Some(primary) = fields.next() else {
            continue;
        };
        let Some(secondary) = fields.next() else {
            return Err(GuideError::InvalidChrMapLine {
                line_no: line_no + 1,
                reason: "expected 2 columns: primary_chrom secondary_chrom".to_string(),
            });
        };
        if fields.next().is_some() {
            return Err(GuideError::InvalidChrMapLine {
                line_no: line_no + 1,
                reason: "expected exactly 2 columns".to_string(),
            });
        }

        if chrmap
            .insert(primary.to_string(), secondary.to_string())
            .is_some()
        {
            return Err(GuideError::InvalidChrMapLine {
                line_no: line_no + 1,
                reason: format!("duplicate primary_chrom: {primary}"),
            });
        }
    }

    Ok(chrmap)
}

#[derive(Debug)]
struct ParsedBedRecord {
    chrom: String,
    strand: ISOMSTRAND,
    interval: GuideInterval,
}

#[derive(Debug)]
pub enum GuideError {
    Io { path: PathBuf, source: io::Error },
    InvalidBedLine { line_no: usize, reason: String },
    InvalidChrMapLine { line_no: usize, reason: String },
}

impl std::fmt::Display for GuideError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GuideError::Io { path, source } => {
                write!(f, "I/O error when loading {}: {source}", path.display())
            }
            GuideError::InvalidBedLine { line_no, reason } => {
                write!(f, "invalid BED line {line_no}: {reason}")
            }
            GuideError::InvalidChrMapLine { line_no, reason } => {
                write!(f, "invalid chrmap line {line_no}: {reason}")
            }
        }
    }
}

impl std::error::Error for GuideError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            GuideError::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}

fn parse_bed_record(line: &str, line_no: usize) -> Result<ParsedBedRecord, GuideError> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 6 {
        return Err(GuideError::InvalidBedLine {
            line_no,
            reason: format!("expected at least 6 columns, got {}", fields.len()),
        });
    }

    let chrom = fields[0].to_string();
    let bed_start = parse_u32_field(fields[1], line_no, "start")?;
    let bed_end = parse_u32_field(fields[2], line_no, "end")?;
    if bed_end <= bed_start {
        return Err(GuideError::InvalidBedLine {
            line_no,
            reason: format!(
                "BED end must be greater than start for half-open interval, got start={} end={}",
                bed_start, bed_end
            ),
        });
    }

    let score = parse_f32_field(fields[4], line_no, "score")?;
    let strand = parse_strand_field(fields[5], line_no)?;

    // BED is 0-based half-open [start, end); convert to 1-based closed [start+1, end].
    let interval = GuideInterval {
        start: bed_start + 1,
        end: bed_end,
        score,
    };

    Ok(ParsedBedRecord {
        chrom,
        strand,
        interval,
    })
}

fn parse_u32_field(raw: &str, line_no: usize, field_name: &str) -> Result<u32, GuideError> {
    raw.parse::<u32>().map_err(|_| GuideError::InvalidBedLine {
        line_no,
        reason: format!("invalid {field_name}: {raw}"),
    })
}

fn parse_f32_field(raw: &str, line_no: usize, field_name: &str) -> Result<f32, GuideError> {
    raw.parse::<f32>().map_err(|_| GuideError::InvalidBedLine {
        line_no,
        reason: format!("invalid {field_name}: {raw}"),
    })
}

fn parse_strand_field(raw: &str, line_no: usize) -> Result<ISOMSTRAND, GuideError> {
    match raw {
        "+" => Ok(ISOMSTRAND::Plus),
        "-" => Ok(ISOMSTRAND::Minus),
        "." => Ok(ISOMSTRAND::Unknown),
        _ => Err(GuideError::InvalidBedLine {
            line_no,
            reason: format!("invalid strand: {raw}"),
        }),
    }
}
