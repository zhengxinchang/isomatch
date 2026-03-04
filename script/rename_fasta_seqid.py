#!/usr/bin/env python3

import argparse
import sys
from typing import TextIO


def rewrite_fasta_seqid(
    input_handle: TextIO,
    output_handle: TextIO,
    sep: str = "|",
) -> None:
    for line in input_handle:
        if not line.startswith(">"):
            output_handle.write(line)
            continue

        line_ending = ""
        if line.endswith("\r\n"):
            header = line[1:-2]
            line_ending = "\r\n"
        elif line.endswith("\n"):
            header = line[1:-1]
            line_ending = "\n"
        else:
            header = line[1:]

        parts = header.split(None, 1)
        seqid = parts[0]
        rest = parts[1] if len(parts) > 1 else ""

        new_seqid = seqid.split(sep, 1)[0] if sep else seqid
        if rest:
            output_handle.write(f">{new_seqid} {rest}{line_ending}")
        else:
            output_handle.write(f">{new_seqid}{line_ending}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rewrite FASTA seqid to the first token split by a separator."
    )
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument(
        "-o",
        "--output",
        help="Output FASTA file, default: stdout",
    )
    parser.add_argument(
        "-s",
        "--sep",
        default="|",
        help='Separator used to split seqid, default: "|"',
    )
    args = parser.parse_args()

    if args.output:
        with open(args.input_fasta, "r", encoding="utf-8") as fin, open(
            args.output, "w", encoding="utf-8"
        ) as fout:
            rewrite_fasta_seqid(fin, fout, args.sep)
    else:
        with open(args.input_fasta, "r", encoding="utf-8") as fin:
            rewrite_fasta_seqid(fin, sys.stdout, args.sep)


if __name__ == "__main__":
    main()
