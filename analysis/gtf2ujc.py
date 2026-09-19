#!/usr/bin/env python3
"""
read a gtf file and convert it to a table file

read all multiple exon transcript and make ujc id based on

chr name, strand, sorted intron chain, 用: 分割

output a table

column 1  transccript_id
column 2  exon number
column 3  ujc

output another table

用ujc 合并多个相同UJC的转录本，输出一个表格，包含
column 1  ujc
column 2  exon number
column 3  transcript number
column 4  transcript id list 逗号分隔

同时输出一个统计信息，包含多少个转录本输入，多少个是多外显子转录本，多少个是单外显子转录本，多少个是多外显子转录本的唯一ujc id

usage: gtf2ujc.py <in.gtf> <out_prefix>
outputs: <out_prefix>.tx_ujc.tsv, <out_prefix>.ujc.tsv, <out_prefix>.stats.tsv
UJC format: chr:strand:intron1_start-intron1_end:intron2_start-intron2_end... (1-based, inclusive)
"""
import re
import sys
from collections import defaultdict

TX_ID_RE = re.compile(r'transcript_id "([^"]+)"')


def read_exons(gtf):
    # transcript_id -> [chrom, strand, [(start, end), ...]]
    txs = {}
    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            tx_id = TX_ID_RE.search(fields[8]).group(1)
            tx = txs.setdefault(tx_id, [fields[0], fields[6], []])
            tx[2].append((int(fields[3]), int(fields[4])))
    return txs


def make_ujc(chrom, strand, exons):
    exons = sorted(exons)
    introns = [f"{e1 + 1}-{s2 - 1}" for (_, e1), (s2, _) in zip(exons, exons[1:])]
    return ":".join([chrom, strand] + introns)


def main(gtf, prefix):
    txs = read_exons(gtf)
    ujc_txs = defaultdict(list)
    n_single = 0

    with open(f"{prefix}.tx_ujc.tsv", "w") as out:
        out.write("transcript_id\texon_number\tujc\n")
        for tx_id, (chrom, strand, exons) in txs.items():
            if len(exons) < 2:
                n_single += 1
                continue
            ujc = make_ujc(chrom, strand, exons)
            ujc_txs[ujc].append(tx_id)
            out.write(f"{tx_id}\t{len(exons)}\t{ujc}\n")

    with open(f"{prefix}.ujc.tsv", "w") as out:
        out.write("ujc\texon_number\ttranscript_number\ttranscript_ids\n")
        for ujc, ids in ujc_txs.items():
            n_exon = len(txs[ids[0]][2])
            out.write(f"{ujc}\t{n_exon}\t{len(ids)}\t{','.join(ids)}\n")

    stats = [
        ("input_transcripts", len(txs)),
        ("multi_exon_transcripts", len(txs) - n_single),
        ("single_exon_transcripts", n_single),
        ("multi_exon_unique_ujc", len(ujc_txs)),
    ]
    with open(f"{prefix}.stats.tsv", "w") as out:
        for key, value in stats:
            out.write(f"{key}\t{value}\n")
            print(f"{key}\t{value}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
