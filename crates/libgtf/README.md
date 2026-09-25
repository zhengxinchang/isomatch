# libgtf

`libgtf` is the library behind isomatch's GTF processing and index format. It
can parse and group GTF transcripts, access indexed FASTA files, build `.isomx`
and `.isoms` files, and read those files without depending on the isomatch CLI.

## Add The Crate

In this workspace:

```toml
[dependencies]
libgtf = { path = "crates/libgtf" }
```

The main APIs are available from these modules:

| Module | Purpose |
| --- | --- |
| crate root | Build an `.isomx`/`.isoms` index pair |
| `libgtf::gtf` | Parse GTF records and stream grouped transcripts |
| `libgtf::fasta` | Read regions from an indexed FASTA file |
| `libgtf::index` | Read index headers, chromosomes, transcripts, and attributes |
| `libgtf::error` | Structured library errors and the shared `Result` alias |

## Build An Index

`build_index` is the simplest entry point. Every supplied FASTA must have an
adjacent `.fai` index.

```rust,no_run
use std::path::PathBuf;

use libgtf::{BuildConfig, build_index};

fn main() -> libgtf::error::Result<()> {
    let report = build_index(&BuildConfig {
        gtf_path: PathBuf::from("sample.gtf.gz"),
        reference_fasta: PathBuf::from("reference.fa"),
        transcript_fasta: None,
        isomx_output: PathBuf::from("sample.isomx"),
        isoms_output: PathBuf::from("sample.isoms"),
        skip_missing_reference_seqids: false,
        temp_dir: None,
    })?;

    println!("indexed {} transcripts", report.transcript_count);
    Ok(())
}
```

### `BuildConfig`

| Field | Meaning |
| --- | --- |
| `gtf_path` | Input GTF path. Plain-text and gzip-compressed input are supported. |
| `reference_fasta` | Reference genome FASTA used to calculate junction and sequence data. |
| `transcript_fasta` | Optional transcript FASTA used to calculate transcript sequence hashes. |
| `isomx_output` | Destination for the transcript index. |
| `isoms_output` | Destination for the compressed original GTF attributes. |
| `skip_missing_reference_seqids` | Skip GTF seqids absent from the reference instead of returning an error. |
| `temp_dir` | Optional parent directory for temporary GTF sorting files. |

`BuildReport` contains transcript and gene totals, skipped records, missing
seqids, strand totals, exon totals, and canonical/non-canonical junction
statistics.

### Build Events

Use `build_index_with_events` when an application needs progress messages. The
callback is synchronous and is called from the build operation. `BuildEvent` is
non-exhaustive, so consumers must include a fallback match arm.

```rust,no_run
use libgtf::{BuildConfig, BuildEvent, build_index_with_events};

fn build(config: &BuildConfig) -> libgtf::error::Result<()> {
    let report = build_index_with_events(config, |event| match event {
        BuildEvent::ProcessingChromosome { name } => {
            eprintln!("processing {name}");
        }
        BuildEvent::MissingReferenceSequence { seqid } => {
            eprintln!("skipping missing reference sequence {seqid}");
        }
        _ => {}
    })?;

    println!("indexed {} genes", report.gene_count);
    Ok(())
}
```

`build_index` performs the same operation with a no-op event callback.

## Parse GTF Data

### Individual Records

`parse_record` borrows fields directly from one GTF line. `attribute` finds an
exact attribute key and supports the quoted and unquoted forms accepted by the
original isomatch parser.

```rust
use libgtf::gtf::{Strand, attribute, parse_record};

let line = "chr1\tsource\texon\t12\t34\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\";";
let record = parse_record(line)?;

assert_eq!(record.seqid, "chr1");
assert_eq!(record.feature, "exon");
assert_eq!(record.strand, Strand::Minus);
assert_eq!(attribute(record.attributes, "transcript_id"), Some("T1"));
# Ok::<(), libgtf::error::Error>(())
```

`GtfRecord` exposes `seqid`, `feature`, `start`, `end`, `strand`, and the raw
`attributes` column. Coordinates are 1-based and inclusive.

### Grouped Transcripts

`GtfReader` groups transcript and exon records, sorts them deterministically,
and returns one `Transcript` at a time. It uses temporary files so that the
whole GTF does not have to be retained in memory.

```rust,no_run
use libgtf::gtf::GtfReader;

fn read_gtf() -> libgtf::error::Result<()> {
    let mut reader = GtfReader::new("sample.gtf.gz")?;

    println!("{} chromosomes", reader.profile().chrom_names.len());
    while let Some(tx) = reader.next()? {
        println!("{}\t{}:{}-{}", tx.tx_id, tx.chrom_id, tx.start, tx.end);
    }
    Ok(())
}
```

Use `GtfReader::new_in(input, temp_parent)` to choose the parent directory for
temporary files. `GtfReader` deliberately does not implement `Iterator` because
reading can fail; call its fallible `next() -> Result<Option<Transcript>>`
method instead.

Useful GTF types and methods:

| API | Purpose |
| --- | --- |
| `GtfProfile` | Chromosome names and IDs, input hash, and input file size |
| `Transcript` | IDs, chromosome, strand, bounds, exons, and original attributes |
| `GtfReader::chrom_name` | Resolve a `ChromID` to its sequence name |
| `GtfReader::transcript_count_excluding` | Count transcripts after excluding chromosome IDs |
| `Strand` | `Plus`, `Minus`, or `Unknown`; converts to/from its index bit value |

## Read FASTA Data

`FastaReader` requires a FASTA index at `<fasta-path>.fai`. Unlike GTF and
index coordinates, FASTA fetch ranges are **0-based and half-open**.

```rust,no_run
use libgtf::fasta::{FaType, FastaReader};

fn fetch_sequence() -> Result<(), libgtf::error::FastaError> {
    let mut fasta = FastaReader::open("reference.fa", FaType::Ref)?;
    let bases = fasta.fetch("chr1", 99, 109, false)?;

    assert_eq!(bases.len(), 10);
    println!("reference contains chr1: {}", fasta.contains("chr1"));
    Ok(())
}
```

`fetch_all` returns an entire sequence. `contains`, `seq_len`, and `seqids`
inspect the FASTA index without fetching bases. `FaType::Ref` and `FaType::Seq`
identify reference and transcript FASTA call sites. The current `trim_chr`
argument to `fetch` and `fetch_all` is reserved for compatibility and does not
alter sequence names.

## Read An Index

`IndexReader` reads `.isomx`; `AttrIndexReader` reads the matching `.isoms`.
Both readers return errors rather than terminating the process.

```rust,no_run
use std::fs::File;

use libgtf::index::{AttrIndexReader, IndexReader};

fn read_index() -> libgtf::error::Result<()> {
    let mut index = IndexReader::open(File::open("sample.isomx")?, 0)?;
    let mut attrs = AttrIndexReader::open("sample.isoms")?;

    for chrom_name in index.chromosome_names().to_vec() {
        let mut chrom = index.get_chromosome_reader(&chrom_name)?;

        while let Some(tx) = chrom.next_record()? {
            let tx_id = tx.source_tx_id(chrom.string_pool());
            let junctions = tx.junction_slice(chrom.junction_pool())?;
            let raw_attributes = attrs.get_attr(tx.tx_idx())?;

            println!(
                "{chrom_name}\t{tx_id}\t{}-{}\t{} junction coordinates\tattributes={}",
                tx.start(),
                tx.end(),
                junctions.len(),
                raw_attributes.is_some(),
            );
        }
    }
    Ok(())
}
```

### Reader API

| API | Purpose |
| --- | --- |
| `IndexReader::load_header` | Decode only the `.isomx` header |
| `IndexReader::open` | Open an `.isomx`; `file_id` is retained on chromosome readers |
| `version`, `md5`, `transcript_count` | Read index metadata |
| `chromosome_names`, `contains_chromosome` | Inspect indexed chromosomes |
| `missing_seqids` | Return GTF seqids skipped while the index was built |
| `get_chromosome_reader` | Open one chromosome block |
| `get_chromosome_readers_map` | Open all chromosome blocks keyed by name |
| `build_txid_index` | Map source transcript IDs to transcript indices |
| `ChromBlockReader::next_record` | Fallibly read the next `TxBase` |
| `ChromBlockReader::reset` | Rewind a chromosome reader to its first transcript |
| `AttrIndexReader::get_attr` | Decompress original attributes by transcript index |

`TxBase` provides transcript metadata through `tx_idx`, `chrom_id`, `start`,
`end`, `strand`, `flags`, `n_exons`, `seq_hash`, `ref_hash`, and
`tx_boundary`. Pool-backed data is resolved with `junction_slice`, `junctions`,
`splice_sites`, `source_tx_id`, and `source_gene_id` using the corresponding
pool from `ChromBlockReader`.

The lower-level `JunctionPool`, `SpliceSitePool`, `StringPool`, their span
types, `TxBaseFlags`, and `TxBoundary` are also public for consumers that need
to inspect or construct the compact index data model. Most applications only
need the readers above.

## Errors And Logging

Public operations return either `libgtf::error::Result<T>` or a specific error:

| Error | Scope |
| --- | --- |
| `Error` | GTF processing, index building, and index reading |
| `FastaError` | FASTA opening, lookup, bounds, and I/O failures |
| `IndexDataError` | Invalid compact values, spans, pools, or encodings |

Errors are structured and retain their underlying sources where applicable.
The library never calls `process::exit` and never initializes a logger. It may
emit diagnostics through the `log` facade; a binary such as isomatch decides
whether and how to install a logger.

Index readers are intended for indexes produced by a compatible libgtf or
isomatch version. Validation of every malformed offset and pool span is not yet
complete, so applications should not treat untrusted index files as fully
validated input.

## Format Compatibility

Current on-disk format versions are exported as:

- `libgtf::index::ISOMX_VERSION` for `.isomx`
- `libgtf::index::ISOMS_VERSION` for `.isoms`
- `libgtf::index::ISOM_GTF_SCHEMA` for the GTF-derived schema

A format change must update the corresponding version. Parsed transcript and
index coordinates remain 1-based and inclusive. Chromosomes are ordered
lexicographically; transcripts are then ordered by chromosome, start, end,
strand, transcript ID, and gene ID.
