# isomatch  <img src="./img/logo.png" align="right" alt="" width=180 />

**Why isomatch?**

Merge

- Supports configurable wobble for splice sites, TSS, and TES, making it well suited for long-read RNA-seq data.
- Considers both junction structure and TSS/TES during transcript merging, rather than relying only on intron chains.
- Uses splice-junction sequence context (GT-AG/GC-AG/AT-AC) to prioritize transcripts during merging.
- Supports guide-based representative TSS/TES selection, enabling integration of high-confidence resources such as refTSS and PolyASite.
- Fully tracks source transcripts after merging.

Classification

- Provides comprehensive transcript classification in seconds to minutes for millions of transcripts.
- Supports SQANTI3-compatible structural classification categories.

Scalability

- Fast and memory-efficient.
- Designed for population-scale analysis across thousands of samples (~10GB memory).

Easy to use

- No dependencies.
- Download and run.

# Subcommands and workflow

- `isomatch index`: build `.isomx` and `.isoms` indexes for a GTF with a reference FASTA. GTF records can be unordered; transcripts must have `transcript` and `exon` records with matching seqids.

- `isomatch merge`: merge transcripts from multiple inputs, auto-rebuilding missing/corrupted/outdated indexes when `--ref-fa` is supplied, and select representative transcripts based on optional third-party evidence.

- `isomatch classify`: classify query transcripts against a reference annotation, auto-rebuilding missing/corrupted/outdated query/reference indexes when `--ref-fa` is supplied.

- `isomatch tools chop`: remove selected GTF attributes, commonly used to strip isomatch-added annotations from output GTF files.

- `isomatch tools valtable`: extract a per-transcript attribute （e.g., TPM）value from source GTFs and assemble it into a matrix aligned to the merged GTF transcript order.

Typical workflow:
1. create indexes with `isomatch index --ref-fa ref.fa`, or let `merge`/`classify` auto-create them.
2. merge transcripts with `isomatch merge --ref-fa ref.fa`, optionally using TSS/TES guide evidence.
3. classify query or merged transcripts with `isomatch classify --ref-fa ref.fa --ref-gtf ref.gtf.gz`.
4. use `isomatch tools` to manipulate GTF outputs, e.g., 

    1. `isomatch tools chop` to remove isomatch-added attributes, 
    2. `isomatch tools valtable` to extract expression values from source GTFs.


# Installation

1. Download the latest release binary from the [GitHub releases page](https://github.com/zhengxinchang/isomatch/releases)
2. Unpack the archive
3. Grant executable permission if needed: Linux: `chmod +x isomatch`, MacOS: `xattr -d com.apple.quarantine isomatch`
4. Run `./isomatch --help` to see usage instructions

# Examples

### Merge
```

# download the isoamtch binary from the latest release page

# optional pre-indexing; merge can auto-index the same inputs if needed
isomatch index --ref-fa ref.fa sample1.gtf.gz
isomatch index --ref-fa ref.fa sample2.gtf.gz
isomatch index --ref-fa ref.fa sample3.gtf.gz

# merge with default parameters (-o takes a prefix, not a filename)
# merge checks <input>.isomx and <input>.isoms; missing/stale indexes are rebuilt
isomatch merge --ref-fa ref.fa -o merged sample1.gtf.gz sample2.gtf.gz sample3.gtf.gz
# outputs: merged.merged.gtf.gz  merged.track.tsv.gz  merged.merged_info.json

# merge with guide-based terminal selection
isomatch merge --ref-fa ref.fa -o merged \
    --guide-tss human.grch38.tss.bed --guide-tes human.grch38.tes.bed \
    sample1.gtf.gz sample2.gtf.gz sample3.gtf.gz

# merge with wobble splice junction matching
isomatch merge --ref-fa ref.fa -o merged \
    -d 3 -a 3 -u 3 \
    -D 5 -A 5 -U 5 \
    sample1.gtf.gz sample2.gtf.gz sample3.gtf.gz
```

### Classify
```
# classify query transcripts against a reference annotation
# classify also checks/rebuilds query and reference indexes
isomatch classify --ref-fa ref.fa --ref-gtf reference.gtf.gz \
    -o query_vs_ref query.gtf.gz
# outputs: query_vs_ref.classification.txt.gz  query_vs_ref.annotated.gtf.gz  query_vs_ref.classify_info.json
```

# How isomatch merge transcripts

## Design Principles

The core idea of isomatch merge is: **the intron chain defines transcript identity; TSS and TES distinguish isoforms sharing the same chain.**

Unlike intron chain collapse only tools, isomatch:

1. Uses splice junction matching (with configurable wobble tolerance) as the primary merge criterion;
2. Further splits transcripts with identical splice chains by TSS/TES distance thresholds;
3. Treats canonical (GT-AG/GC-AG/AT-AC) and non-canonical splice sites separately, with independent wobble parameters;
4. Supports third-party terminal evidence (e.g., refTSS, PolyA site databases) to guide representative TSS/TES selection.

---

<details>
<summary><strong>Read details</strong></summary>


## Pipeline Overview

```
Multiple GTF inputs
        │
        ▼
[1] isomatch index: each GTF + reference FASTA → .isomx and .isoms indexes
        │
        ▼
[2] K-way merge + genomic coordinate overlap → Super Cluster (per-chromosome locus)
        │
        ▼
[3] Group by (strand, exon count) → Junction Cluster
        │
        ├─ canonical (ALLC)
        ├─ non-canonical (PRTC / NOTC)
        └─ mono-exon (MONO)
        │
        ▼
[4] Splice junction wobble matching → initial groups
        │
        ▼
[5] Terminal Refine: further split groups by TSS/TES distance
        │
        ▼
[6] Representative transcript selection (splice junction / TSS / TES policy)
        │
        ▼
Output GTF (with ISOM_SRC / ISOM_COUNT / ISOM_REPR_POLICY annotations)
```

---

## Stage-by-Stage Details

### Stage 1: Index Building (isomatch index)

Each input GTF is preprocessed together with the reference FASTA into `.isomx` and `.isoms` index files. GTF input may be unordered; isomatch aggregates exon records by `transcript_id` and sorts transcripts internally. Each transcript is stored with its chromosome, strand, coordinates, exon count, splice junctions, and transcript type:

- `ALLC`: all splice sites are canonical (GT-AG, GC-AG, AT-AC)
- `PRTC`: partially canonical splice sites
- `NOTC`: no canonical splice sites
- `MONO`: single-exon transcript, no splice sites

---

### Stage 2: Super Cluster Construction

For each chromosome, transcripts from all input files are streamed in coordinate order and grouped into **super clusters** — sets of transcripts that overlap each other on the genome, equivalent to a locus.

---

### Stage 3: Junction Cluster Grouping

Within each super cluster, transcripts are further partitioned by `(strand, n_exons)` into **junction clusters**.

This guarantees merge homogeneity: transcripts with different exon counts or on different strands are never merged together.

---

### Stage 4: Splice Junction Wobble Matching (core merge step)

#### 4a. Canonical transcript merging

Transcripts of type `ALLC` are grouped using a **Union-Find** algorithm with wobble matching:

- For any pair of transcripts, each corresponding junction's donor and acceptor coordinates are compared;
- If all junctions fall within the wobble thresholds, the two transcripts are unioned into the same group.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--wob-d` | Max donor-site deviation (canonical) | 0 bp |
| `--wob-a` | Max acceptor-site deviation (canonical) | 0 bp |
| `--wob-u` | Splice-site wobble for unstranded transcripts (canonical) | 3 bp |

By default, canonical transcripts require **exact junction matches** (donor/acceptor wobble = 0).

#### 4b. Non-canonical transcript merging

`PRTC` / `NOTC` transcripts first attempt to be **absorbed into an existing canonical group**: each non-canonical transcript is compared against every canonical group's representative junction chain using the nc-wobble parameters, and assigned to the best-matching group (lowest total junction deviation).

Non-canonical transcripts that cannot be absorbed are then grouped among themselves using the same Union-Find approach with nc-wobble parameters:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--wob-d-nc` | Max donor-site deviation (non-canonical) | 3 bp |
| `--wob-a-nc` | Max acceptor-site deviation (non-canonical) | 3 bp |
| `--wob-u-nc` | Splice-site wobble for unstranded transcripts (non-canonical) | 3 bp |

#### 4c. Mono-exon transcript merging

`MONO` transcripts are merged based on **reciprocal overlap**:

$$\text{reciprocal overlap} = \min\!\left(\frac{|A \cap B|}{|A|},\ \frac{|A \cap B|}{|B|}\right)$$

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--mono-ovlp` | Minimum reciprocal overlap threshold | 0.9 |

---

### Stage 5: Terminal Refinement

After initial junction-based grouping, groups are further **split by TSS and TES distance**. This ensures that transcripts sharing the same intron chain but with substantially different transcription start or end sites are kept as distinct isoforms.

Algorithm: transcripts within a group are sorted by TSS/TES; the first transcript sets the anchor. Each subsequent transcript is compared to the current anchor — if the distance exceeds the threshold, a new group begins:

$$
\lvert \mathrm{TSS}_{\mathrm{curr}} - \mathrm{TSS}_{\mathrm{anchor}} \rvert \leq \tau_{\mathrm{TSS}}
\quad \text{and/or, depending on mode} \quad
\lvert \mathrm{TES}_{\mathrm{curr}} - \mathrm{TES}_{\mathrm{anchor}} \rvert \leq \tau_{\mathrm{TES}}
$$

Here, $\tau_{\mathrm{TSS}}$ and $\tau_{\mathrm{TES}}$ correspond to `--tss-wob` and `--tes-wob` (or `--tss-wob-nc` and `--tes-wob-nc` for non-canonical transcripts).

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--terminal-refine` | Refine mode (`none`/`tss`/`tes`/`both`) for canonical transcripts | `both` |
| `--tss-wob` | Max TSS deviation within a group (canonical) | 50 bp |
| `--tes-wob` | Max TES deviation within a group (canonical) | 50 bp |
| `--terminal-refine-nc` | Refine mode for non-canonical transcripts | `both` |
| `--tss-wob-nc` | Max TSS deviation within a group (non-canonical) | 50 bp |
| `--tes-wob-nc` | Max TES deviation within a group (non-canonical) | 50 bp |

---

### Stage 6: Representative Transcript Selection

Each merged group outputs one representative transcript. Selection is performed independently for three dimensions:

#### Splice junctions

Each intron's representative coordinates are chosen independently:

| Policy | Description |
|--------|-------------|
| `major` (default) | Most frequent junction pair; falls back to `longer` on tie |
| `longer` | Shortest intron (longest flanking exon span) |
| `shorter` | Longest intron (shortest flanking exon span) |

Note: for splice junctions, `longer` means **longer exons** (shorter intron), consistent with the "longer transcript" convention.

#### TSS and TES

| Policy | TSS meaning | TES meaning |
|--------|-------------|-------------|
| `major` (default) | Most frequent TSS; falls back to most upstream on tie | Most frequent TES; falls back to most downstream on tie |
| `longer` | Most upstream TSS | Most downstream TES |
| `shorter` | Most downstream TSS | Most upstream TES |

Controlled by `--tss-policy` and `--tes-policy` (both default to `major`).

#### Guide-based terminal selection

When `--guide-tss` and/or `--guide-tes` BED files are provided (e.g., refTSS, PolyA site databases), terminal selection is biased toward positions supported by external evidence:

1. Each TSS/TES position in the group is queried against the guide database within `±flank` bp (default 10 bp);
2. Candidates with the highest number of guide evidence hits are retained;
3. Majority vote is applied among those candidates; ties fall back to the `longer` policy.

When a guide hit is used, the corresponding policy column records one of three values:

| Value | Meaning |
|-------|---------|
| `guide_definitive` | All guide-supported candidates pointed to the same position; guide alone determines the choice |
| `guide_majority` | Multiple guide-supported candidates; majority vote picked one winner |
| `guide_longer` | Guide-supported candidates were tied after majority vote; fell back to the `longer` rule |

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--guide-tss` | TSS evidence BED file | none |
| `--guide-tes` | TES evidence BED file | none |
| `--guide-tss-flank` | TSS evidence query flank radius | 10 bp |
| `--guide-tes-flank` | TES evidence query flank radius | 10 bp |

#### Mono-exon boundary selection

Mono-exon transcripts are first grouped by reciprocal overlap using `--mono-ovlp`, then their representative boundaries are chosen with the same guide-aware `--tss-policy` / `--tes-policy` terminal selection described above.

---

</details>

## Output Files

For a run with `-o <prefix>`, three files are written:

| File | Description |
|------|-------------|
| `<prefix>.merged.gtf.gz` | Merged transcripts in GTF format (gzip-compressed) |
| `<prefix>.track.tsv.gz` | One-to-one mapping of merged → source transcripts (gzip-compressed TSV) |
| `<prefix>.merged_info.json` | Run statistics (source file count, merged transcript counts, guide usage, etc.) |

In `<prefix>.merged_info.json`, `*_pct` values are percentages on a 0-100 scale rounded to four decimal places.

---

## Output GTF Format

Each merged transcript record in the output GTF includes the following extra attributes on the `transcript` line:

| Attribute | Description |
|-----------|-------------|
| `ISOM_EXONS` | Number of exons |
| `ISOM_COUNT` | Number of source transcripts merged into this record |
| `ISOM_SRC` | `\|`-separated list of source transcripts, each formatted as `S{sample_index}:{tx_id}:{start}:{end}:{tx_type}:{donor_diff}:{acceptor_diff}:{exon_diffs}` |
| `ISOM_REPR_POLICY` | Representative selection policies as `SJ_POLICY:TSS_POLICY:TES_POLICY`; `SJ_POLICY` is `NA` for mono-exon transcripts |

The `exon_diffs` field in `ISOM_SRC` records only exons that differ from the representative, in the format `(exon_number,left_offset,right_offset)`. Exons with no difference are omitted and shown as `no_diff`.

Policy values for all policy fields: `major`, `longer`, `shorter`, `guide_definitive`, `guide_majority`, `guide_longer`.

---

## Track TSV Format

`<prefix>.track.tsv.gz` has one row per (merged transcript, source transcript) pair. Columns:

| Column | Description |
|--------|-------------|
| `merged_tx_id` | Merged transcript ID (`ISOMT_{n}`) |
| `merged_gene_id` | Merged gene ID (`ISOMG_{n}`) |
| `merged_start` | Merged transcript left genomic boundary |
| `merged_end` | Merged transcript right genomic boundary |
| `merged_strand` | Strand (`+`, `-`, or `.`) |
| `merged_exon_num` | Number of exons in the merged transcript |
| `junction_policy` | Splice junction representative policy used; `NA` for mono-exon transcripts |
| `tss_policy` | TSS representative policy used (strand-aware) |
| `tes_policy` | TES representative policy used (strand-aware) |
| `src_tx_count_in_merged_group` | Total number of source transcripts merged into this record |
| `src_tx_id` | Source transcript ID |
| `src_gene_id` | Source gene ID |
| `total_donor_diff` | Sum of donor-site deviations between source and representative junctions |
| `total_acceptor_diff` | Sum of acceptor-site deviations between source and representative junctions |
| `exon_diff` | Per-exon coordinate differences in `(exon_number,left_offset,right_offset)` format; `no_diff` if identical |

Note: `junction_policy`, `tss_policy`, and `tes_policy` are the same for all rows sharing the same `merged_tx_id`, and `ISOM_REPR_POLICY` now uses the same strand-aware `SJ_POLICY:TSS_POLICY:TES_POLICY` ordering.

# How isomatch classify transcripts

## Design Principles

The core idea of isomatch classify is: **compare each query transcript against the reference annotation, choose the best transcript-level hit, then report SQANTI3-style structural categories with sequence context.**

isomatch classify:

1. Uses reference junctions, splice sites, exons, and gene spans as catalog evidence;
2. Classifies multi-exon transcripts primarily by ordered splice-junction relationships;
3. Handles mono-exon transcripts by exon overlap, containment, and intron-retention evidence;
4. Adds strand-aware TSS/TES differences, optional CAGE/polyA guide evidence, and FASTA-derived TES sequence context.

---

<details>
<summary><strong>Read details</strong></summary>

## Classify Pipeline Overview

```
Query GTF + Reference GTF + Reference FASTA
        │
        ▼
[1] Load query indexes (.isomx/.isoms), reference index (.isomx), and reference genome sequence
        │
        ▼
[2] For each query transcript: collect global reference catalog evidence
    (known junctions, known splice sites, known introns, antisense overlaps)
        │
        ▼
[3] Find overlapping reference transcripts on the same strand
        │
        ▼
[4] Compare query vs each reference transcript
    ├─ exact junction chain      → FSM
    ├─ consecutive subset chain  → ISM
    ├─ known junction/site usage → NIC or NNC
    └─ exon/gene overlap         → genic classes
        │
        ▼
[5] Select the primary hit per gene, then the best overall hit
        │
        ▼
[6] Refine multi-gene, antisense, genic-intron, and intergenic cases
        │
        ▼
[7] Add FASTA sequence context around TES and optional CAGE/polyA guide evidence
        │
        ▼
Output classification table, annotated GTF, and summary JSON
```

---

## Classification Rules

For multi-exon query transcripts, isomatch compares the ordered junction chain to each overlapping reference transcript:

| Relationship | Main category |
|--------------|---------------|
| Same junction chain | `full-splice_match` |
| Query junction chain is a consecutive subset of reference | `incomplete-splice_match` |
| Novel combination using known junctions or splice sites from the same gene | `novel_in_catalog` |
| At least one novel splice site | `novel_not_in_catalog` |
| Exon overlap without splice-chain evidence | `genic` |

FSM subcategories are determined by strand-aware TSS/TES differences. The match threshold is controlled by `--fsm-end-match-bp` (default 50 bp).

For mono-exon query transcripts, isomatch uses exon overlap against reference transcripts:

| Relationship | Main category |
|--------------|---------------|
| Overlaps a mono-exon reference transcript | `full-splice_match` |
| Fully contained in a reference exon | `incomplete-splice_match` |
| Spans a known reference intron | `novel_in_catalog` |
| Otherwise overlaps a same-strand gene | `genic` |

If no same-strand hit is found, isomatch falls back to `antisense`, `genic_intron`, or `intergenic`, based on opposite-strand overlap and known intron containment.

---

</details>

## Classify Output

For a run with `-o <prefix>`, the following files are written:

| File | Description |
|------|-------------|
| `<prefix>.classification.txt.gz` | Per-transcript classification and sequence-context metrics |
| `<prefix>.annotated.gtf.gz` | Query transcript/exon GTF annotated with classification attributes |
| `<prefix>.classify_info.json` | Run statistics (query/reference counts, category/subcategory counts, guide support, etc.) |

The annotated GTF preserves query transcript/exon structures and adds or refreshes `ISOM_REF_TX_ID`, `ISOM_REF_GENE_ID`, `ISOM_REF_GENE_NAME`, `ISOM_CATEGORY`, and `ISOM_SUBCATEGORY`.

Key output fields include query transcript information, structural category/subcategory, selected reference gene/transcript, TSS/TES differences, matched junction/exon counts, canonical splice status, CAGE/polyA guide support, and downstream TES A content.

In `<prefix>.classify_info.json`, `*_pct` values and `structural_category_pct` values are percentages on a 0-100 scale rounded to four decimal places.


# isomatch Tools

## isomatch tools chop

`chop` removes attributes from a GTF while preserving the first eight GTF columns and comments. By default, it removes `ISOM_*` attributes and keeps standard identifiers such as `gene_id` and `transcript_id`.

```
isomatch tools chop -o cleaned merged.merged.gtf.gz
# outputs: cleaned.chopped.gtf.gz
```

| Item | Description |
|------|-------------|
| Input | GTF or gzip-compressed GTF |
| Output | `<prefix>.chopped.gtf.gz` |

Key parameters:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-o, --out` | Output prefix | required |
| `-m, --mode` | Attribute removal mode: `isomatch` removes only `ISOM_*` attributes; `all` removes all attributes except kept ones | `isomatch` |
| `-k, --keep` | Comma-separated extra attributes to keep | `gene_id,transcript_id` are always kept |
| `-c, --keep-check-case` | Match attribute names case-sensitively | case-insensitive matching |

## isomatch tools valtable

`valtable` reads a merged GTF produced by `isomatch merge`, looks up the source transcripts recorded in each `ISOM_SRC` attribute, extracts a specified GTF attribute value from the original source GTFs, and writes a matrix with one row per merged transcript and one column per source file.

```
isomatch tools valtable \
    -m merged.merged.gtf.gz \
    -o expr \
    -a TPM \
    sample1.gtf.gz sample2.gtf.gz sample3.gtf.gz
# outputs: expr.valtable.tsv.gz  expr.valtable_stats.json
```

Each source GTF is matched to its corresponding sample in the merged GTF by filename. Source files whose filenames do not appear in the merged GTF `##ISOM <SAMPLE>` headers are skipped with a warning. Transcripts in the source GTF that are not recorded in the merged GTF are counted in the stats but do not affect the output.

| Item | Description |
|------|-------------|
| Input | Merged GTF (gzip or plain) + one or more source GTFs (gzip or plain) |
| Output | `<prefix>.valtable.tsv.gz`, `<prefix>.valtable_stats.json` |

Key parameters:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-m, --merged-gtf` | Merged GTF produced by `isomatch merge` | required |
| `-o, --out` | Output prefix | required |
| `-a, --attr` | GTF attribute name to extract from source transcripts (case-sensitive) | required |
| `-d, --default-val` | Value used when a transcript has no entry in a source file | `0.0` |
| (positional) | One or more source GTF files | required |

`<prefix>.valtable.tsv.gz` is a tab-separated matrix: the first column is `transcript_id` (merged transcript IDs in their original order), and each subsequent column is named after the source GTF filename.

`<prefix>.valtable_stats.json` records the number of merged transcripts, the list of source samples found in the merged GTF header, and per-file extraction counts (`extracted`, `attr_missing`, `src_tx_not_in_merged`).