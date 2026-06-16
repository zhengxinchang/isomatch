use crate::core::tx_strand::ISOMSTRAND;
use crate::index::gtf::parse_gtf_attr_value;
use crate::utils::open_file_bufread;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

// bed format for guide tss tes
// chromosome      start   end     ID       score   strand
// chr1    16013   16020   rfhg_1.1        1       -

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RegionType {
    Tss,
    Tes,
    Gene,
}

pub type ChromMap = FxHashMap<String, String>;

#[derive(Debug, Clone, PartialEq)]
pub struct MyRegion {
    pub start: u32, // 1-based closed
    pub end: u32,   // 1-based closed
    pub score: f32,
    pub id: String,
    pub name: String,
    pub strand: ISOMSTRAND,
}

impl MyRegion {
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
    intervals: Vec<MyRegion>,
    max_len: u32,
}

impl ChromGuideIndex {
    pub fn intervals(&self) -> &[MyRegion] {
        &self.intervals
    }

    pub fn max_len(&self) -> u32 {
        self.max_len
    }

    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }

    pub fn query_overlaps(&self, pos: u32) -> Vec<&MyRegion> {
        self.query_overlaps_with_flank(pos, 0)
    }

    pub fn query_overlaps_with_flank(&self, pos: u32, flank: u32) -> Vec<&MyRegion> {
        let start = pos.saturating_sub(flank);
        let end = pos.saturating_add(flank);
        self.query_overlaps_range(start, end)
    }

    fn query_overlaps_range(&self, start: u32, end: u32) -> Vec<&MyRegion> {
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
pub struct RegionDb {
    guide_type: RegionType,
    // bed_chroms: HashSet<String>,
    by_chrom_strand: FxHashMap<(String, ISOMSTRAND), ChromGuideIndex>,
    chrmap: Option<ChromMap>,
}

impl RegionDb {
    pub fn from_bed_path<P: AsRef<Path>>(
        path: P,
        guide_type: RegionType,
        chrmap_path: &Option<P>,
    ) -> Result<Self, RegionError> {
        let path = path.as_ref();
        let file = File::open(path).map_err(|err| RegionError::Io {
            path: path.to_path_buf(),
            source: err,
        })?;
        let reader = BufReader::new(file);
        // Self::from_bed_reader(reader, guide_type)
        // let mut bed_chroms = HashSet::default();

        let mut grouped: FxHashMap<(String, ISOMSTRAND), Vec<MyRegion>> = FxHashMap::default();

        for (line_no, line_result) in reader.lines().enumerate() {
            let raw_line = line_result.map_err(|err| RegionError::Io {
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
                let max_len = intervals.iter().map(MyRegion::len).max().unwrap_or(0);
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

    pub fn from_gtf_gene<P: AsRef<Path>>(
        path: P,
        guide_type: RegionType,
        chrmap_path: &Option<P>,
    ) -> Result<Self, RegionError> {
        let path = path.as_ref();
        let mut reader = open_file_bufread(path).map_err(|err| RegionError::Io {
            path: path.to_path_buf(),
            source: err,
        })?;
        let mut grouped: FxHashMap<(String, ISOMSTRAND), Vec<MyRegion>> = FxHashMap::default();
        let mut line = String::new();
        let mut line_no = 0usize;

        while reader.read_line(&mut line).map_err(|err| RegionError::Io {
            path: path.to_path_buf(),
            source: err,
        })? != 0
        {
            line_no += 1;
            let raw = line.trim_end();
            if raw.is_empty() || raw.starts_with('#') {
                line.clear();
                continue;
            }

            let fields: Vec<&str> = raw.split('\t').collect();
            if fields.len() < 9 {
                return Err(RegionError::InvalidGtfLine {
                    line_no,
                    reason: format!("expected 9 columns, got {}", fields.len()),
                });
            }
            if fields[2] != "gene" {
                line.clear();
                continue;
            }

            let start = parse_gtf_u32_field(fields[3], line_no, "start")?;
            let end = parse_gtf_u32_field(fields[4], line_no, "end")?;
            let strand = parse_gtf_strand_field(fields[6], line_no)?;
            let id = parse_gtf_attr_value(fields[8], "gene_id").unwrap_or_default();
            let name = parse_gtf_attr_value(fields[8], "gene_name").unwrap_or_else(|| id.clone());

            grouped
                .entry((fields[0].to_string(), strand))
                .or_default()
                .push(MyRegion {
                    start,
                    end,
                    score: 0.0,
                    id,
                    name,
                    strand,
                });
            line.clear();
        }

        let by_chrom_strand = grouped
            .into_iter()
            .map(|(key, mut intervals)| {
                intervals.sort_by_key(|interval| interval.start);
                let max_len = intervals.iter().map(MyRegion::len).max().unwrap_or(0);
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
            chrmap,
        })
    }

    pub fn guide_type(&self) -> RegionType {
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

    pub fn query_overlaps(&self, chrom: &str, strand: ISOMSTRAND, pos: u32) -> Vec<&MyRegion> {
        self.get_index(chrom, strand)
            .map(|index| index.query_overlaps(pos))
            .unwrap_or_default()
    }

    pub fn query_overlaps_range_all_strands(
        &self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Vec<&MyRegion> {
        [ISOMSTRAND::Plus, ISOMSTRAND::Minus, ISOMSTRAND::Unknown]
            .into_iter()
            .flat_map(|strand| {
                self.get_index(chrom, strand)
                    .map(|index| index.query_overlaps_range(start, end))
                    .unwrap_or_default()
            })
            .collect()
    }

    pub fn query_overlaps_with_flank(
        &self,
        chrom: &str,
        strand: &ISOMSTRAND,
        pos: u32,
        flank: u32,
    ) -> Vec<&MyRegion> {
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

pub fn load_chrmap_path<P: AsRef<Path>>(path: P) -> Result<ChromMap, RegionError> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|err| RegionError::Io {
        path: path.to_path_buf(),
        source: err,
    })?;
    let reader = BufReader::new(file);
    let mut chrmap = FxHashMap::default();

    for (line_no, line_result) in reader.lines().enumerate() {
        let raw_line = line_result.map_err(|err| RegionError::Io {
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
            return Err(RegionError::InvalidChrMapLine {
                line_no: line_no + 1,
                reason: "expected 2 columns: primary_chrom secondary_chrom".to_string(),
            });
        };
        if fields.next().is_some() {
            return Err(RegionError::InvalidChrMapLine {
                line_no: line_no + 1,
                reason: "expected exactly 2 columns".to_string(),
            });
        }

        if chrmap
            .insert(primary.to_string(), secondary.to_string())
            .is_some()
        {
            return Err(RegionError::InvalidChrMapLine {
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
    interval: MyRegion,
}

#[derive(Debug)]
pub enum RegionError {
    Io { path: PathBuf, source: io::Error },
    InvalidBedLine { line_no: usize, reason: String },
    InvalidGtfLine { line_no: usize, reason: String },
    InvalidChrMapLine { line_no: usize, reason: String },
}

impl std::fmt::Display for RegionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RegionError::Io { path, source } => {
                write!(f, "I/O error when loading {}: {source}", path.display())
            }
            RegionError::InvalidBedLine { line_no, reason } => {
                write!(f, "invalid BED line {line_no}: {reason}")
            }
            RegionError::InvalidGtfLine { line_no, reason } => {
                write!(f, "invalid GTF line {line_no}: {reason}")
            }
            RegionError::InvalidChrMapLine { line_no, reason } => {
                write!(f, "invalid chrmap line {line_no}: {reason}")
            }
        }
    }
}

impl std::error::Error for RegionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            RegionError::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}

fn parse_bed_record(line: &str, line_no: usize) -> Result<ParsedBedRecord, RegionError> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 6 {
        return Err(RegionError::InvalidBedLine {
            line_no,
            reason: format!("expected at least 6 columns, got {}", fields.len()),
        });
    }

    let chrom = fields[0].to_string();
    let bed_start = parse_u32_field(fields[1], line_no, "start")?;
    let bed_end = parse_u32_field(fields[2], line_no, "end")?;
    if bed_end <= bed_start {
        return Err(RegionError::InvalidBedLine {
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
    let interval = MyRegion {
        start: bed_start + 1,
        end: bed_end,
        score,
        id: fields[3].to_string(),
        name: fields[3].to_string(),
        strand,
    };

    Ok(ParsedBedRecord {
        chrom,
        strand,
        interval,
    })
}

fn parse_gtf_u32_field(raw: &str, line_no: usize, field_name: &str) -> Result<u32, RegionError> {
    raw.parse::<u32>().map_err(|_| RegionError::InvalidGtfLine {
        line_no,
        reason: format!("invalid {field_name}: {raw}"),
    })
}

fn parse_gtf_strand_field(raw: &str, line_no: usize) -> Result<ISOMSTRAND, RegionError> {
    match raw {
        "+" => Ok(ISOMSTRAND::Plus),
        "-" => Ok(ISOMSTRAND::Minus),
        "." => Ok(ISOMSTRAND::Unknown),
        _ => Err(RegionError::InvalidGtfLine {
            line_no,
            reason: format!("invalid strand: {raw}"),
        }),
    }
}

fn parse_u32_field(raw: &str, line_no: usize, field_name: &str) -> Result<u32, RegionError> {
    raw.parse::<u32>().map_err(|_| RegionError::InvalidBedLine {
        line_no,
        reason: format!("invalid {field_name}: {raw}"),
    })
}

fn parse_f32_field(raw: &str, line_no: usize, field_name: &str) -> Result<f32, RegionError> {
    raw.parse::<f32>().map_err(|_| RegionError::InvalidBedLine {
        line_no,
        reason: format!("invalid {field_name}: {raw}"),
    })
}

fn parse_strand_field(raw: &str, line_no: usize) -> Result<ISOMSTRAND, RegionError> {
    match raw {
        "+" => Ok(ISOMSTRAND::Plus),
        "-" => Ok(ISOMSTRAND::Minus),
        "." => Ok(ISOMSTRAND::Unknown),
        _ => Err(RegionError::InvalidBedLine {
            line_no,
            reason: format!("invalid strand: {raw}"),
        }),
    }
}
