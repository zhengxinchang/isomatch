# isomatch

Isomatch: evidence-baesd transcirpts merging and classification.

evidence based merge
    1. not solely based on coordinates, 
    2. consider the canonical splice site pattern,
    3. transcirpt sequence for further comparison, 
    4. TES/TSS motifs from FANTOM5 or CAGE-seq data, 
    5. population evidence (isopedia derived evidence)

multi classification system anntation:
    1. GFFcompare classification
    2. SQANTI3 classification
    3. IsoSeq classification


other feature
    1. full level sample tracking
    2. large-scale processing capability (thousands of samples in GTF format)



# How isomatch merge transcripts

1. treat the splice junction chain as first priority evidence.
1.1 all canonical transcripts will have higher score than non-all canonical tx
1.2 for a junction in all transcripts that sightly differnt, most common of the splice junction will be choosen.
1.3 non-all canonical tx will be attached into the canonical transcirpt
1.4 the rest of non-all canonical tx will be further merged based on the splice junction wobble and tss tes threshold 

## multi exon & plus minus strand

    wobble sj

    wobble tss tes 

    splice site canonical or not 

    guided splice junction
    guided tss tes

    ref squence hash for rescue tx 

    sample sequence hash for ?? 

## multi exon unstrand

## mono exon plus mimnus

## mono exon unstrand


guided 会默认比其他policy 更优先，
如果guidesd失败，则还是会选择policy，所以policy是一直存在的。


case1
合并long read RNA，都是canonical 有1-2bp的差异，则合并？ 