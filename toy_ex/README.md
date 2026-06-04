# isomatch toy example

- `merge` with 3 input GTF files
- `classify` of the merged output against a reference annotation
- automatic indexing, canonical / non-canonical junction handling, terminal shifts, mono-exon matching, antisense, and intergenic cases

## Files

| File | Description |
| --- | --- |
| `ref.fa.gz` | A 1000 bp toy chromosome named `chrToy`, compressed with BGZF/gzip |
| `ref.fa.gz.fai`, `ref.fa.gz.gzi` | Random-access indexes created by `samtools faidx` for the compressed FASTA |
| `ref.gtf.gz` | Reference annotation with multi-exon GeneA/GeneB transcripts and one mono-exon GeneC transcript |
| `sample1.gtf.gz` | Merge input 1, containing FSM, ISM, mono-exon, and antisense examples |
| `sample2.gtf.gz` | Merge input 2, containing terminal shift, NIC, NNC, and mono-exon-contained examples |
| `sample3.gtf.gz` | Merge input 3, containing a 1 bp non-canonical wobble case, far TSS, mono-exon overlap, intron-retention, and intergenic examples |

The GTF inputs are gzip/BGZF-compressed. The commands below automatically create `.isomx` and `.isoms` indexes when needed.

## Merge

Run from the repository root:

```bash
cargo build --release
rm -rf toy_ex/out
mkdir -p toy_ex/out

./target/release/isomatch merge \
  --ref-fa toy_ex/ref.fa.gz \
  -o toy_ex/out/toy_merge \
  toy_ex/sample1.gtf.gz toy_ex/sample2.gtf.gz toy_ex/sample3.gtf.gz
```

Expected output files:

```text
toy_ex/out/toy_merge.merged.gtf.gz
toy_ex/out/toy_merge.track.tsv.gz
toy_ex/out/toy_merge.merged_info.json
```

Expected merge summary:

```json
{
  "source_files": 3,
  "total_tx_cnt": 13,
  "merged_tx_cnt": 10,
  "merged_multi_exons_tx_cnt": 5,
  "merged_mono_exon_tx_cnt": 5,
  "tss_guide_cnt": 0,
  "tes_guide_cnt": 0,
  "merged_tx_by_source_cnt": {
    "1": 7,
    "2": 3
  }
}
```

## Classify

Classify the merged output:

```bash
./target/debug/isomatch classify \
  --ref-fa toy_ex/ref.fa.gz \
  --ref-gtf toy_ex/ref.gtf.gz \
  -o toy_ex/out/toy_classify \
  toy_ex/out/toy_merge.merged.gtf.gz
```

Expected output files:

```text
toy_ex/out/toy_classify.classification.txt.gz
toy_ex/out/toy_classify.annotated.gtf.gz
toy_ex/out/toy_classify.classify_info.json
```

Expected classify category counts:

```json
{
  "antisense": 1,
  "full-splice_match": 3,
  "incomplete-splice_match": 2,
  "intergenic": 1,
  "novel_in_catalog": 2,
  "novel_not_in_catalog": 1
}
```

## Clean Up

```bash
rm -rf toy_ex/out
rm -f toy_ex/*.gtf.gz.isomx toy_ex/*.gtf.gz.isoms toy_ex/*.gtf.gz.isomx.info.json
```
