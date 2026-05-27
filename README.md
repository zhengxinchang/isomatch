# isomatch




Isomatch: evidence-baesd transcirpts merging and classification.

evidence based merge

1. not only consider intron chain, but also consider the TSS/TES
2. configurable splice junction wobble and TSS/TES threshold for merging,
3. select representative transcripts based on third party evidence such as refTSS and PolyAsites

multi classification system anntation:

1. provide SQANTI3 style classification code.
2. for merged gtf from isomatch merge, not only report difference between reference and merged transcripts, but also report difference between the refernce transcripts and the original transcripts.

large-scale processing capability:
1. support thousands of GTF files

# subcomands and workflow

isomatch index: build index for GTF files together with a reference FASTA, required by the rest of the subcommands, and can be reused for multiple runs of merge and classify.

isomatch merge: merge transcripts from multiple indexed inputs, and select representative transcripts based on third party evidence(optional)

isomatch classify: classify transcripts in a GTF file based on reference annotation, and report the classification code for each transcript.

Typical workflow:
1. create index for all input GTF files using `isomatch index --ref-fa ref.fa`
2. merge transcripts from multiple indexed inputs using `isomatch merge`, and select representative transcripts based on third party evidence(optional)
3. classify merged transcripts using `isomatch classify`, and report the classification code for each transcript.


# Examples

### Merge
```
# build index, GTF must be sorted at chromosome level.
# use bedtools sort -i input.gtf > sorted_input.gtf to sort if needed.
isomatch index --ref-fa ref.fa sample1.gtf.gz
isomatch index --ref-fa ref.fa sample2.gtf.gz
isomatch index --ref-fa ref.fa sample3.gtf.gz

# merge with default parameters (-o takes a prefix, not a filename)
# merge uses the same input paths that were indexed and reads <input>.isomx
isomatch merge -o merged  sample1.gtf.gz sample2.gtf.gz sample3.gtf.gz
# outputs: merged.merged.gtf.gz  merged.track.tsv.gz  merged.merged_info.json

# merge with guide-based terminal selection
isomatch merge -o merged \
    --guide-tss human.guide.tss.bed --guide-tes human.guide.tes.bed \
    sample1.gtf.gz sample2.gtf.gz sample3.gtf.gz

# merge with wobble splice junction matching
isomatch merge -o merged \
    -d 3 -a 3 -u 3 \
    -D 5 -A 5 -U 5 \
    sample1.gtf.gz sample2.gtf.gz sample3.gtf.gz
```

# How isomatch merge transcripts

## Design Principles

The core idea of isomatch merge is: **the intron chain defines transcript identity; TSS and TES distinguish isoforms sharing the same chain.**

Unlike simple coordinate-overlap collapse tools, isomatch:

1. Uses splice junction matching (with configurable wobble tolerance) as the primary merge criterion;
2. Further splits transcripts with identical splice chains by TSS/TES distance thresholds;
3. Treats canonical (GT-AG/GC-AG/AT-AC) and non-canonical splice sites separately, with independent wobble parameters;
4. Supports third-party terminal evidence (e.g., refTSS, PolyA site databases) to guide representative TSS/TES selection.

---

## Pipeline Overview

```
Multiple GTF inputs
        │
        ▼
[1] isomatch index: each GTF + reference FASTA → .isomx binary index
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

Each input GTF is preprocessed together with the reference FASTA into a `.isomx` binary index file. Each transcript is stored with its chromosome, strand, coordinates, exon count, splice junctions, and transcript type:

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

## Output Files

For a run with `-o <prefix>`, three files are written:

| File | Description |
|------|-------------|
| `<prefix>.merged.gtf.gz` | Merged transcripts in GTF format (gzip-compressed) |
| `<prefix>.track.tsv.gz` | One-to-one mapping of merged → source transcripts (gzip-compressed TSV) |
| `<prefix>.merged_info.json` | Run statistics (source file count, merged transcript counts, guide usage, etc.) |

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
