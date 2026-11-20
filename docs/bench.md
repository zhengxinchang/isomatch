# isomatch bench

`isomatch bench` evaluates one or more query GTFs against a reference GTF at the
transcript level and reports TP, FP, FN, precision, recall and F1 for each query.

```
isomatch bench --ref-fa ref.fa -g reference.gtf -o bench_out query1.gtf query2.gtf
```

## Presets

| | `--preset default` | `--preset junction` |
|---|---|---|
| Evaluated records | transcripts of each input GTF | each input collapsed on its own by exact junction chain (mono-exon: exact start and end) |
| Multi-exon candidate | junctions within `w`, TSS/TES within `--tss-wob`/`--tes-wob` | junctions within `w`; ends ignored |
| Mono-exon candidate | TSS/TES thresholds and reciprocal overlap | reciprocal overlap only |

- `default` asks whether the full-length transcript is right; `junction` asks
  whether the splicing structure is right, like the gffcompare intron chain
  level. Partial chains never match in either preset.
- Short-read assemblies have reliable junctions but unreliable ends: use
  `junction`. Long-read GTFs: `junction` for structure, `default` for full
  length.
- `--tss-wob`/`--tes-wob` only apply to `default`; giving them with `junction`
  is an error.
- Collapsing runs `isomatch merge` on each input alone with every splice wobble
  0, no terminal refinement and `--mono-ovlp 1.0`. It does not depend on
  `--wob`, so the totals are the same for every `w`.

## Sample labels

Outputs refer to the reference as `REF` and to queries as `Q1`, `Q2`, ... in
command-line order. Each query is evaluated independently against the full
reference.

The intermediate merge step numbers inputs `S1` (reference), `S2`, ... as usual
for `isomatch merge`; those IDs only appear inside `initial_group/`
(`S<k+1>` is `Q<k>`).

## Matching

- Query and reference transcripts are compared only within the same chromosome,
  strand and exon count, and their spans must overlap.
- `default`: every candidate needs TSS and TES within `--tss-wob` / `--tes-wob`
  (50 bp unless given).
- Multi-exon: the largest junction coordinate difference must be within `w`.
  Results are reported for every `w = 0..=--wob`; `w = 0` (identical junction
  chain) is the primary result.
- Mono-exon: reciprocal overlap of at least `--mono-ovlp`; results do not depend
  on `w`.
- Pairs are one-to-one, accepted greedily in the order (junction difference,
  TSS difference, TES difference, query ID, reference ID). With `junction`,
  multi-exon pairs skip the TSS/TES differences; mono-exon pairs keep them.
- Totals are the transcripts of the evaluated GTF accepted by `isomatch index`
  (`<gtf>.index_info.json`): the input GTF, or its collapsed GTF with `junction`.

## Comparison with gffcompare

gffcompare is the usual tool for transcript-level accuracy. The two tools define
a match differently, so their numbers for the same inputs can differ.

| | `isomatch bench` | gffcompare (transcript level) |
|---|---|---|
| Pairing | one-to-one; every query and reference transcript is counted once, as TP or FP/FN | per query class codes; several queries can get `=` for the same reference transcript in `.tmap`, while `.stats` counts matching transcripts once per reference |
| Redundant transcripts | kept on both sides (`default`); collapsed on both sides by exact structure (`junction`) | duplicate reference transcripts are collapsed; query duplicates are kept unless `-D` is given |
| Denominators | all transcripts of the evaluated query and reference GTFs, stored in the outputs | Query mRNAs / Reference mRNAs in `.stats`, after gffcompare's own filtering |
| Coordinates | used as given | exons separated by very short introns are merged first |
| Ends | `--tss-wob`/`--tes-wob` separately (`default`), or ignored (`junction`) | one threshold, `-e`, applied with `--strict-match` |
| Tolerance | results for every junction wobble `w = 0..=--wob` in one run | exact intron chain |
| Per-transcript output | `tp-query`, `fp`, `tp-ref`, `fn` GTFs that keep the source IDs | `.tmap` / `.refmap` with class codes |
| Inputs | reference FASTA required (for indexing) | GTF only |

What this changes in practice:

- **Redundant transcripts.** A GTF that contains transcripts with identical
  structure, compared with itself, gives F1 = 1 with bench. gffcompare
  collapses the duplicate reference transcripts but keeps the duplicate
  queries, so its transcript-level precision and F1 are below 1 although the
  files are identical. You can check this on any such GTF.
- **Counting only `=` hides false positives.** A common way to score a tool
  from gffcompare output counts the `=` predictions and de-duplicates their
  `cmp_ref` IDs. Predictions with other class codes then never enter the
  precision denominator. bench counts every prediction, so precision can be
  lower, especially for tools that report many partial or redundant models.
- **Very short introns.** Predictions whose junction chain differs from the
  reference only by an intron of a few bp can be `=` in gffcompare, because it
  merges such exons first, but are FP with bench, which compares the
  coordinates as given.
- **Gaps between tools can change.** Because the denominators and match rules
  differ, tools are not affected equally: re-scoring the same predictions with
  bench can widen or narrow the differences between tools, not only shift all
  values together.

These differences come from different definitions and from how gffcompare
output is post-processed, not only from the matching algorithm. Neither F1 is
"the correct one"; compare tools only with numbers computed the same way.

Use bench when you need one-to-one counts, every prediction in the
denominator, matching sensitivity across `w`, or per-transcript TP/FP/FN files.
gffcompare remains the better fit for class codes (`j`, `k`, `c`, ...), exon-,
intron- and locus-level statistics, or when no reference FASTA is available.

## Output layout

```
bench_out/
├── params.json          # all bench parameters
├── results_all.json     # metrics per query and w, with full input summaries
├── results_summary.json # the same metrics with input paths and totals only
├── results_summary.tsv  # one row per query and w: all-record counts and ratios
├── collapsed/           # junction only: one collapsed input per sample
│   ├── REF/
│   │   ├── REF.merged.gtf.gz          # evaluated records
│   │   ├── REF.merged.gtf.gz.isomx, ...  # its index
│   │   └── REF.track.tsv.gz, ...      # collapsed ID -> original IDs
│   └── Q1/ ...
├── initial_group/       # intermediate merge output used for candidate grouping
│   └── bench.merged.gtf.gz, bench.track.tsv.gz, ...
├── Q1/
│   ├── tp-query.gtf.gz  # Q1 transcripts matched to a reference transcript
│   ├── fp.gtf.gz        # Q1 transcripts without a match
│   ├── tp-ref.gtf.gz    # reference transcripts matched by Q1
│   └── fn.gtf.gz        # reference transcripts not matched by Q1
└── Q2/ ...
```

Every bench output directory contains a `.ISOMATCH_BENCH` marker file. Output
directory rules:

| Output directory | without `--overwrite` | with `--overwrite` |
|---|---|---|
| missing | created | created |
| empty | used | used |
| non-empty, has `.ISOMATCH_BENCH` | error | contents deleted, then run |
| non-empty, no `.ISOMATCH_BENCH` | error | error, nothing deleted |

With `--overwrite`, nothing is deleted either if the output directory contains an
input GTF, the reference FASTA or its `.fai`, or the current working directory,
or if it is the root or home directory. The directory itself is kept, and
symlinks inside it are removed, not followed.

## results_all.json

A JSON object with:

- `preset`: `default` or `junction`.
- `reference`: `sample_id` (`REF`), `input_path`, `total`, `gtf_stats` and
  `collapsed`.
- `query_count`: number of query GTFs.
- `queries`: one entry per query in command-line order, with the same fields as
  `reference` plus `results`, a list with one item per `w`.

Each `results` item holds `tp_query`, `tp_ref`, `fp`, `fn`, `precision`,
`recall` and `f1` for all records, and the same seven fields again in
`multi_exon` and `mono_exon` for each class alone. Mono-exon results do not
depend on `w`.

`gtf_stats` is the input file's `<gtf>.index_info.json` written by `isomatch
index` (transcript, strand, exon and splice-site counts, missing seqids).
`collapsed` is `null` with `default`; with `junction` it has the collapsed GTF
`path` and its own `gtf_stats`. `total` is the `transcript_count` of the
evaluated GTF. Ratios have six decimals; a ratio with a zero denominator is
`null`.

## results_summary.json

The same `results` as `results_all.json` (counts, precision, recall and F1 for
all, multi-exon and mono-exon records), without the index summaries:

- `preset`, and `evaluated_on`: `input` (default) or `collapsed` (junction), the
  GTFs the metrics are computed on.
- `params`: `wob`, `tss_wob`, `tes_wob` (`null` with junction) and `mono_ovlp`.
  `params.json` has every parameter.
- `reference`, `query_count` and `queries`. Each input has:

| Field | Meaning |
|---|---|
| `sample_id` | `REF`, `Q1`, ... |
| `input_path` | command-line GTF |
| `evaluated_path` | GTF the metrics are computed on: `input_path` with default, the collapsed GTF with junction |
| `input_total` | transcripts of `input_path` accepted by `isomatch index` |
| `total` | transcripts of `evaluated_path`: the denominator of precision (query) and recall (reference) |
| `multi_exon_total`, `mono_exon_total` | denominators of the `multi_exon` and `mono_exon` results |

Queries also have `results`.

## results_summary.tsv

Tab-separated, one row per query and `w`, for all records (no multi/mono-exon
split): `sample_id`, `input_path`, `w`, `tp_query`, `tp_ref`, `fp`, `fn`,
`precision`, `recall`, `f1`. A ratio with a zero denominator is `NA`. The `w = 0`
rows are also printed to the log at the end of the run.

## Status GTFs

The four files under each `Q#/` use the `w = 0` result. Every transcript of the
query appears in exactly one of `tp-query`/`fp`, and every reference transcript
in exactly one of `tp-ref`/`fn`.

Records are rebuilt from the merge output: each transcript has one `transcript`
line and its `exon` lines, with the `gene_id` and `transcript_id` of the
evaluated GTF,
`ISOM_TX_ID` (merged transcript in `initial_group/bench.merged.gtf.gz`) and
`exon_number` in genomic order. Other original attributes are not kept.
With `junction` the IDs and coordinates are those of the collapsed GTF
(`ISOMT_*`); its `ISOM_SRC` attribute or `track.tsv.gz` gives the original
transcripts.

Header example (`Q2/fp.gtf.gz`):

```
##ISOM <VERSION> version="0.8.0"; program="isomatch"; schema="1"
##ISOM <STATUS_BENCH> id="fp"; side="query"; w="0"; preset="default"; description="query transcripts of Q2 without a one-to-one reference match at w=0";
##ISOM <SAMPLE> id="REF"; input="reference.gtf"; role="reference"; total="602025";
##ISOM <SAMPLE> id="Q2"; input="query2.gtf"; role="query"; total="413280";
##ISOM <FORMAT> ID="ISOM_TX_ID"; Description="merged transcript ID in initial_group/bench.merged.gtf.gz of the group holding this transcript";
##ISOM <COMMAND> cmd=isomatch bench ...
```

- `STATUS_BENCH` names the status (`tp-query`, `fp`, `tp-ref`, `fn`), the side
  (`query` or `reference`), the `w` used and the preset.
- `SAMPLE` lists the reference and this query. In bench status GTFs the `id` is
  a bench label (`REF`, `Q#`), not a merge `S#` ID; `role` and `total` (index
  transcript count of the evaluated GTF) are bench fields. `input` is the
  evaluated GTF; with `junction` it is the collapsed GTF and `raw_input` is the
  input it came from:

  ```
  ##ISOM <SAMPLE> id="Q2"; input="bench_out/collapsed/Q2/Q2.merged.gtf.gz"; raw_input="query2.gtf"; role="query"; total="8406";
  ```
- `COMMAND` may appear more than once: tools that rewrite the file (for example
  `isomatch classify`) append their own line after the existing ones.
- Parameter defaults follow the version in `VERSION`; `params.json` has the full
  values.
