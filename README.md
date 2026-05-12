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


## multi exon unstrand

## mono exon plus mimnus

## mono exon unstrand


guide 模式不是annotate，将merge之后的转录本合并到参考转录本，而是用第三方证据例如fantom5的tss tes motif，或者isopedia的population-scale的splice junction evidence来指导merge的过程。


guided 会默认比其他policy 更优先，
如果guidesd失败，则还是会选择policy，所以policy是一直存在的。


guide bed 的格式

chr start end stand score

case1
合并long read RNA，都是canonical 有1-2bp的差异，则合并？


Merge 和 Annotate 都不应该修改转录本的任何坐标，
guided模式在merge中的作用是帮助选择representatibve sites。

merge 和 refine 是为了确定哪些转录本应该是一个， guide 是为了帮助选择代表性的转录本。这两个是分开的步骤。


smallexon rescue 应该是一个单独的功能或者整合到correct（refine）中，因为涉及到对转录本的修改。

what is missing in merge
stats
guide 
ISOM_COUNT --> ISOM_TX_CNT plus ISOM_SMPLE_CNT

ISOM_EXONS 直接显示数字，而不是擅自分类

gffcompare 在进行群体级别的合并的时候存在内存问题，如果是steam模式，则问题是所有的TSS TES都是来自第一个参考GTF。
isomatch, 可以在多个参数微调的情况下，使用kwaymerge方法从所有的样本中获得TSS TES的证据，避免内存问题，同时也能充分利用所有样本的证据。

isomatch 的merge 不改变原始输入的转录本坐标，而是通过选择代表性的转录本来进行合并。

isomatch refine 的功能是对任何GTF根据外部信息进行校正，例如小外显子救援，或者根据外部证据调整转录本的TSS TES位置。

isomatch annotate 的功能是对转录本进行分类注释，例如根据gffcompare的分类，或者根据sqanti3的分类，或者根据isoseq的分类。

isomatch compare 