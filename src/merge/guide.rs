use crate::core::tx_strand::ISOMSTRAND;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

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
        let file = File::open(path).map_err(GuideError::Io)?;
        let reader = BufReader::new(file);
        // Self::from_bed_reader(reader, guide_type)
        // let mut bed_chroms = HashSet::default();

        let mut grouped: FxHashMap<(String, ISOMSTRAND), Vec<GuideInterval>> = FxHashMap::default();

        for (line_no, line_result) in reader.lines().enumerate() {
            let raw_line = line_result.map_err(GuideError::Io)?;

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
    let file = File::open(path).map_err(GuideError::Io)?;
    let reader = BufReader::new(file);
    let mut chrmap = FxHashMap::default();

    for (line_no, line_result) in reader.lines().enumerate() {
        let raw_line = line_result.map_err(GuideError::Io)?;

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
    Io(io::Error),
    InvalidBedLine { line_no: usize, reason: String },
    InvalidChrMapLine { line_no: usize, reason: String },
}

impl std::fmt::Display for GuideError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GuideError::Io(err) => write!(f, "I/O error: {err}"),
            GuideError::InvalidBedLine { line_no, reason } => {
                write!(f, "invalid BED line {line_no}: {reason}")
            }
            GuideError::InvalidChrMapLine { line_no, reason } => {
                write!(f, "invalid chrmap line {line_no}: {reason}")
            }
        }
    }
}

impl std::error::Error for GuideError {}

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

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn unique_temp_path(prefix: &str, suffix: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("{prefix}-{nanos}.{suffix}"))
    }

    fn write_temp_file(prefix: &str, suffix: &str, contents: &str) -> PathBuf {
        let path = unique_temp_path(prefix, suffix);
        fs::write(&path, contents).unwrap();
        path
    }

    fn write_temp_bed(contents: &str) -> PathBuf {
        write_temp_file("isomatch-guide", "bed", contents)
    }

    fn write_temp_chrmap(contents: &str) -> PathBuf {
        write_temp_file("isomatch-chrmap", "tsv", contents)
    }

    fn sample_bed_contents() -> &'static str {
        concat!(
            "chromosome\tstart\tend\tID\tscore\tstrand\n",
            "chr1\t9\t12\ta\t1\t+\n",
            "chr1\t19\t22\tb\t2\t+\n",
            "chr1\t29\t40\tc\t3\t-\n",
            "chr1\t100\t180\tlong\t4\t+\n",
            "chr2\t14\t18\td\t5\t+\n",
        )
    }

    #[test]
    fn loads_bed_and_converts_coordinates() {
        let bed_path = write_temp_bed(sample_bed_contents());
        let db = GuideDb::from_bed_path(&bed_path, GuideBEDType::Tss, &None::<&PathBuf>).unwrap();
        let plus = db.get_index("chr1", ISOMSTRAND::Plus).unwrap();
        let minus = db.get_index("chr1", ISOMSTRAND::Minus).unwrap();

        let _ = fs::remove_file(&bed_path);

        assert_eq!(db.guide_type(), GuideBEDType::Tss);
        assert_eq!(db.len(), 3);
        assert_eq!(plus.max_len(), 80);
        assert!(plus.intervals().contains(&GuideInterval {
            start: 10,
            end: 12,
            score: 1.0,
        }));
        assert_eq!(
            minus.intervals(),
            &[GuideInterval {
                start: 30,
                end: 40,
                score: 3.0,
            }]
        );
    }

    #[test]
    fn query_overlaps_returns_point_hits() {
        let bed_path = write_temp_bed(sample_bed_contents());
        let db = GuideDb::from_bed_path(&bed_path, GuideBEDType::Tss, &None::<&PathBuf>).unwrap();

        let hits = db.query_overlaps("chr1", ISOMSTRAND::Plus, 11);

        let _ = fs::remove_file(&bed_path);

        assert_eq!(hits.len(), 1);
        assert_eq!(
            *hits[0],
            GuideInterval {
                start: 10,
                end: 12,
                score: 1.0,
            }
        );
        assert!(db.query_overlaps("chr1", ISOMSTRAND::Minus, 11).is_empty());
        assert!(db.query_overlaps("chr2", ISOMSTRAND::Plus, 11).is_empty());
    }

    #[test]
    fn query_overlaps_with_flank_expands_range() {
        let bed_path = write_temp_bed(sample_bed_contents());
        let db = GuideDb::from_bed_path(&bed_path, GuideBEDType::Tss, &None::<&PathBuf>).unwrap();

        let point_hits = db.query_overlaps("chr1", ISOMSTRAND::Plus, 16);
        let flank_hits = db.query_overlaps_with_flank("chr1", &ISOMSTRAND::Plus, 16, 4);

        let _ = fs::remove_file(&bed_path);

        assert!(point_hits.is_empty());
        assert_eq!(flank_hits.len(), 2);
        assert_eq!(flank_hits[0].start, 10);
        assert_eq!(flank_hits[0].end, 12);
        assert_eq!(flank_hits[1].start, 20);
        assert_eq!(flank_hits[1].end, 22);
    }

    #[test]
    fn query_overlaps_with_flank_keeps_long_interval_candidates() {
        let bed_path = write_temp_bed(sample_bed_contents());
        let db = GuideDb::from_bed_path(&bed_path, GuideBEDType::Tss, &None::<&PathBuf>).unwrap();
        let index = db.get_index("chr1", ISOMSTRAND::Plus).unwrap();

        let hits = index.query_overlaps_with_flank(181, 1);

        let _ = fs::remove_file(&bed_path);

        assert_eq!(hits.len(), 1);
        assert_eq!(
            *hits[0],
            GuideInterval {
                start: 101,
                end: 180,
                score: 4.0,
            }
        );
    }

    #[test]
    fn invalid_bed_line_returns_error() {
        let bed_path =
            write_temp_bed("chromosome\tstart\tend\tID\tscore\tstrand\nchr1\tstart\t12\ta\t1\t+\n");
        let err =
            GuideDb::from_bed_path(&bed_path, GuideBEDType::Tes, &None::<&PathBuf>).unwrap_err();
        let _ = fs::remove_file(&bed_path);

        match err {
            GuideError::InvalidBedLine { line_no, .. } => assert_eq!(line_no, 2),
            other => panic!("unexpected error: {other}"),
        }
    }

    #[test]
    fn loads_chrmap_and_always_skips_first_line() {
        let path = write_temp_chrmap("primary_chrom\tsecondary_chrom\n1\tchr1\nMT\tchrM\n");
        let map = load_chrmap_path(&path).unwrap();
        let _ = fs::remove_file(path);
        assert_eq!(map.len(), 2);
        assert_eq!(map.get("1").map(String::as_str), Some("chr1"));
        assert_eq!(map.get("MT").map(String::as_str), Some("chrM"));
    }

    #[test]
    fn query_with_chrmap_uses_secondary_chrom_after_exact_miss() {
        let bed_path = write_temp_bed(sample_bed_contents());
        let path = write_temp_chrmap("primary_chrom\tsecondary_chrom\n1\tchr1\n");
        let db = GuideDb::from_bed_path(&bed_path, GuideBEDType::Tss, &Some(&path)).unwrap();
        let hits = db.query_overlaps_with_flank("1", &ISOMSTRAND::Plus, 16, 4);

        let _ = fs::remove_file(&bed_path);
        let _ = fs::remove_file(path);

        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0].start, 10);
        assert_eq!(hits[1].start, 20);
    }

    #[test]
    fn duplicate_primary_chrom_in_chrmap_returns_error() {
        let path = write_temp_chrmap("primary_chrom\tsecondary_chrom\n1\tchr1\n1\tchr01\n");
        let err = load_chrmap_path(&path).unwrap_err();
        let _ = fs::remove_file(path);

        match err {
            GuideError::InvalidChrMapLine { line_no, .. } => assert_eq!(line_no, 3),
            other => panic!("unexpected error: {other}"),
        }
    }
}
