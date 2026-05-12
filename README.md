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

isomatch 子命令体系（修订版）

index — 对GTF建立.isi二进制索引（B+ tree），供下游子命令调用。也可由其他子命令按需自动触发，但显式暴露便于预处理大文件和复用。
merge — 多样本合并。k-way merge，从所有样本收集TSS/TES证据，选代表性转录本，不改变原始坐标。解决gffcompare在群体级别的内存问题和evidence偏差问题。
refine — 坐标级校正。根据外部证据调整TSS/TES位置、小外显子救援、strand分配。Scope严格限定在坐标修正和stand等基本信息，不涉及功能注释。接受统一的bed格式证据。
classify — 1-vs-ref结构分类。Soft matching，支持configurable wobble tolerance，输出per-transcript classification code（兼容gffcompare class codes和SQANTI3 FSM/ISM/NIC/NNC体系）。服务对象：生物学家和annotation pipeline。
compare — N-vs-N交叉比较。Binary matching（identical structure，yes/no），支持--match-mode参数（strict=全坐标一致，junction=内部剪接位点一致，chain=exon chain topology一致）。输出overlap矩阵（TSV）和UpSet-compatible intersection集合。不做可视化。服务对象：工具开发者做横向比较。
bench — 1-vs-truth准确性评估。Binary matching，match-mode参数同compare。输出TP/FP/FN/precision/recall/F1，支持per-gene和per-transcript两个粒度。输出干净TSV，不做可视化。服务对象：工具开发者和benchmarking consortium（LRGASP等）。

核心架构关系： index是预处理层；classify是核心比对引擎（soft matching）；compare和bench共享一套独立的binary matching引擎；merge和refine是转录组操作层。classify与compare/bench的比对逻辑是分开的——前者输出分类码，后者只输出boolean。



TODO
monoexon guide 模式
简化grouppgir的代码现在充满match

