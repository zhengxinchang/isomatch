
# What you can do with isomatch

Suppose you start with 20 independent RNA-seq GTF files from two biological
conditions, for example 10 condition-A samples and 10 condition-B samples:

```
A1.gtf.gz A2.gtf.gz ... A10.gtf.gz
B1.gtf.gz B2.gtf.gz ... B10.gtf.gz
```

isomatch lets you move from per-sample transcript calls to a unified,
sample-aware isoform catalog, then reuse that catalog for downstream comparison,
annotation, expression, and provenance checks.

## Need: inspect each GTF before merging

Use `isomatch index` to build indexes and write per-file statistics.

```
isomatch index --ref-fa ref.fa A1.gtf.gz
```

The generated `<input>.index_info.json` summarizes the input GTF, including
transcript counts, strand distribution, mono-exon versus multi-exon composition,
canonical versus non-canonical splice-site composition, and missing reference
sequence IDs.

## Need: build one unified isoform catalog from all 20 samples

Use `isomatch merge` on all GTFs. This collapses matching transcripts across
samples while preserving source provenance.

```
isomatch merge --ref-fa ref.fa -o merged \
    A1.gtf.gz A2.gtf.gz ... A10.gtf.gz \
    B1.gtf.gz B2.gtf.gz ... B10.gtf.gz
```

Key outputs:

| Output | Use it for |
|--------|------------|
| `merged.merged.gtf.gz` | Unified transcript catalog with `ISOM_COUNT`, `ISOM_SAMPLE_CNT`, `ISOM_SAMPLE_FREQ`, and source-tracking attributes |
| `merged.track.tsv.gz` | One-to-one mapping from each merged transcript back to its source transcript IDs |
| `merged.present_absent.tsv.gz` | Presence/absence matrix showing which merged isoforms occur in which samples |
| `merged.merged_info.json` | Merge-level summary statistics |

## Need: compare isoform presence between condition A and condition B

Use the `merged.present_absent.tsv.gz` matrix from `isomatch merge`. Each row is
a merged isoform and each sample column records whether that isoform was observed
in the original GTF.

This table is the simplest starting point for questions such as:

- isoforms found only in condition A or only in condition B
- isoforms recurrent in most samples of one condition
- isoforms shared by both conditions but with different sample frequencies

## Need: connect merged isoforms back to the original transcript calls

Use `merged.track.tsv.gz` from `isomatch merge`, or inspect the `ISOM_SRC`
attribute in `merged.merged.gtf.gz`.

This answers provenance questions such as:

- which original transcript IDs were merged into `ISOMT_123`
- which source GTFs contributed to a merged isoform
- how much each source transcript differs from the representative splice
  junctions

## Need: attach TPM, FPKM, CPM, or another GTF attribute to the merged catalog

Use `isomatch tools valtable` with the merged GTF and the original source GTFs.

```
isomatch tools valtable \
    -m merged.merged.gtf.gz \
    -o tpm \
    -a TPM \
    A1.gtf.gz A2.gtf.gz ... A10.gtf.gz \
    B1.gtf.gz B2.gtf.gz ... B10.gtf.gz
```

The output `tpm.valtable.tsv.gz` is a merged-transcript-by-sample matrix aligned
to the unified catalog. Change `-a TPM` to `FPKM`, `CPM`, or any other
transcript-level attribute present in your source GTFs.

## Need: identify known, novel, antisense, or intergenic isoforms

Use `isomatch classify` to compare the unified catalog against a reference
annotation.

```
isomatch classify --ref-fa ref.fa --ref-gtf reference.gtf.gz \
    -o merged_vs_ref merged.merged.gtf.gz
```

Key outputs:

| Output | Use it for |
|--------|------------|
| `merged_vs_ref.classification.txt.gz` | SQANTI3-style structural category table |
| `merged_vs_ref.annotated.gtf.gz` | Unified GTF annotated with reference gene/transcript IDs and classification labels |
| `merged_vs_ref.classify_info.json` | Category counts and classification summary statistics |

## Need: mark overlapping reference genes without full transcript classification

Use `isomatch tools mark` when you only need gene-overlap labels on transcript
records.

```
isomatch tools mark \
    --ref-gtf reference.gtf.gz \
    -o merged_marked \
    merged.merged.gtf.gz
```

The output `merged_marked.mark.gtf.gz` adds `ISOM_OVLP_GENE` to transcript
records, and `merged_marked.mark_dup_gene.tsv.gz` reports transcripts that
overlap multiple reference genes.

## Need: analysis the overlap between source GTFs

use the `merged.present_absent.tsv.gz` matrix to plot upset plots, Venn diagrams, or other visualizations of isoform sharing across samples.

