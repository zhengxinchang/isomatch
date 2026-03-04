# isomatch

Isomatch compares, merges, and annotates transcripts using sequence-context comparison, motif-aware correction, and population-scale processing.

Key features:

- **Sequence-context comparison**: Compares transcripts by their actual sequence, resolving misalignment artifacts caused by short exons.
- **Motif-aware correction**: Corrects splice sites using canonical motif information, improving the accuracy of isoform structure inference.
- **Long-read native**: Leverages full-length transcript sequences from long-read RNA-seq for more faithful comparison.
- **Dual classification system**: Supports both GFFcompare and SQANTI3 classification schemes in a single run.
- **Population-scale processing**: Handles datasets spanning thousands of samples in GTF format.



