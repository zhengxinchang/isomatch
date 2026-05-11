

# TSS evidence bed

http://reftss.clst.riken.jp/datafiles/current/human/refTSS_v4.1_human_coordinate.hg38.bed.txt.gz
```
cut -f 1-6 atlas.clusters.2.0.GRCh38.96.bed |sed 's/^/chr/' >guide.tes.bed
```

# TES evidence bed
https://polyasite.unibas.ch/download/atlas/2.0/GRCm38.96/atlas.clusters.2.0.GRCm38.96.bed.gz
```
cut -f 1-6 atlas.clusters.2.0.GRCh38.96.bed > guide.tes.bed
```


