#!/usr/bin/env python3

import argparse
import gzip
import json
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


def open_text_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_gtf_attributes(attrs: str) -> tuple[str, str]:
    tx_id = ""
    gene_id = ""

    for attr in attrs.split(";"):
        attr = attr.strip()
        if not attr:
            continue
        if attr.startswith("transcript_id"):
            tx_id = extract_attr_value(attr)
        elif attr.startswith("gene_id"):
            gene_id = extract_attr_value(attr)
        if tx_id and gene_id:
            break

    return tx_id, gene_id


def extract_attr_value(attr: str) -> str:
    q_start = attr.find('"')
    if q_start != -1:
        q_end = attr.find('"', q_start + 1)
        if q_end != -1:
            return attr[q_start + 1 : q_end]
    fields = attr.split()
    return fields[1] if len(fields) > 1 else ""


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ATCGNatcgn", "TAGCNtagcn")
    return seq.translate(table)[::-1].upper()


def normalized_site(site: str, strand: int) -> str:
    if strand == 1:
        return reverse_complement(site)
    return site.upper()


def splice_site_pack(left: str, right: str, strand: int) -> int:
    left = normalized_site(left, strand)
    right = normalized_site(right, strand)

    def code(site: str) -> int:
        if site == "GT":
            return 0
        if site == "AG":
            return 1
        if site == "GC":
            return 2
        if site == "AT":
            return 3
        if site == "AC":
            return 4
        return 5

    left_code = code(left)
    right_code = code(right)
    if strand == 1:
        return (right_code << 4) | left_code
    return (left_code << 4) | right_code


def is_canonical_pair(left: str, right: str, strand: int) -> bool:
    packed = splice_site_pack(left, right, strand)
    donor = packed >> 4
    acceptor = packed & 0x0F
    return (
        (donor == 0 and acceptor == 1)
        or (donor == 2 and acceptor == 1)
        or (donor == 3 and acceptor == 4)
    )


@dataclass(frozen=True)
class FaiEntry:
    length: int
    offset: int
    line_bases: int
    line_width: int


class IndexedFasta:
    def __init__(self, fasta_path: Path):
        self.fasta_path = fasta_path
        self.fai_path = Path(str(fasta_path) + ".fai")
        self.entries = self._load_fai()
        self.handle = fasta_path.open("rb")

    def _load_fai(self) -> dict[str, FaiEntry]:
        entries = {}
        with self.fai_path.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                name, length, offset, line_bases, line_width, *_ = raw_line.rstrip("\n").split("\t")
                entries[name] = FaiEntry(
                    length=int(length),
                    offset=int(offset),
                    line_bases=int(line_bases),
                    line_width=int(line_width),
                )
        return entries

    def fetch(self, chrom: str, start_0: int, end_1: int) -> str:
        entry = self.entries.get(chrom)
        if entry is None:
            raise KeyError(f"Chromosome {chrom!r} not found in FASTA index")
        if not (0 <= start_0 <= end_1 <= entry.length):
            raise ValueError(
                f"Requested interval [{start_0}, {end_1}) outside chromosome {chrom} length {entry.length}"
            )

        pos = start_0
        chunks = []
        while pos < end_1:
            line_idx = pos // entry.line_bases
            line_pos = pos % entry.line_bases
            take = min(end_1 - pos, entry.line_bases - line_pos)
            file_offset = entry.offset + line_idx * entry.line_width + line_pos
            self.handle.seek(file_offset)
            chunks.append(self.handle.read(take).decode("ascii"))
            pos += take
        return "".join(chunks).upper()

    def close(self) -> None:
        self.handle.close()


@dataclass
class TranscriptRecord:
    chrom: str
    strand: int
    tx_id: str
    gene_id: str
    exons: list[tuple[int, int]]

    def add_exon(self, start: int, end: int) -> None:
        self.exons.append((start, end))

    def sorted_exons(self) -> list[tuple[int, int]]:
        return sorted(self.exons, key=lambda exon: exon[0])


@dataclass
class IndexStats:
    transcript_count: int = 0
    gene_count: int = 0
    plus_strand_count: int = 0
    minus_strand_count: int = 0
    mono_exon_count: int = 0
    multi_exon_count: int = 0
    junction_count: int = 0
    canonical_junction_count: int = 0
    noncanonical_junction_count: int = 0
    canonical_junction_ratio: float = 0.0


def iter_transcripts(gtf_path: Path) -> Iterable[TranscriptRecord]:
    current = None

    with open_text_maybe_gzip(gtf_path) as handle:
        for raw_line in handle:
            if not raw_line or raw_line.startswith("#"):
                continue

            parts = raw_line.rstrip("\n").split("\t")
            if len(parts) < 9:
                raise ValueError(f"Invalid GTF line with fewer than 9 columns: {raw_line.rstrip()}")

            feature = parts[2]
            if feature != "exon":
                continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            if start > end:
                continue

            strand = 1 if parts[6] == "-" else 0
            tx_id, gene_id = parse_gtf_attributes(parts[8])

            if current is None:
                current = TranscriptRecord(chrom=chrom, strand=strand, tx_id=tx_id, gene_id=gene_id, exons=[])
                current.add_exon(start, end)
                continue

            if tx_id != current.tx_id:
                yield current
                current = TranscriptRecord(chrom=chrom, strand=strand, tx_id=tx_id, gene_id=gene_id, exons=[])

            current.add_exon(start, end)

    if current is not None:
        yield current


def count_canonical_junctions(tx: TranscriptRecord, fasta: IndexedFasta) -> int:
    exons = tx.sorted_exons()
    if len(exons) <= 1:
        return 0

    tx_start_0 = exons[0][0] - 1
    tx_end_1 = exons[-1][1]
    reference_seq = fasta.fetch(tx.chrom, tx_start_0, tx_end_1)

    exon_offsets = [(start - exons[0][0], end - exons[0][0] + 1) for start, end in exons]
    canonical = 0

    for left_exon, right_exon in zip(exon_offsets, exon_offsets[1:]):
        lstart = left_exon[1]
        lend = left_exon[1] + 2
        rstart = right_exon[0] - 2
        rend = right_exon[0]
        left = reference_seq[lstart:lend]
        right = reference_seq[rstart:rend]
        if is_canonical_pair(left, right, tx.strand):
            canonical += 1

    return canonical


def compute_expected_stats(gtf_path: Path, reffa_path: Path) -> IndexStats:
    stats = IndexStats()
    gene_ids = set()
    fasta = IndexedFasta(reffa_path)

    try:
        for tx in iter_transcripts(gtf_path):
            exons = tx.sorted_exons()
            exon_count = len(exons)
            canonical_junction_count = count_canonical_junctions(tx, fasta)

            stats.transcript_count += 1
            gene_ids.add(tx.gene_id)
            if tx.strand == 1:
                stats.minus_strand_count += 1
            else:
                stats.plus_strand_count += 1

            if exon_count <= 1:
                stats.mono_exon_count += 1
                continue

            junction_count = exon_count - 1
            stats.multi_exon_count += 1
            stats.junction_count += junction_count
            stats.canonical_junction_count += canonical_junction_count
            stats.noncanonical_junction_count += junction_count - canonical_junction_count
    finally:
        fasta.close()

    stats.gene_count = len(gene_ids)
    if stats.junction_count:
        stats.canonical_junction_ratio = stats.canonical_junction_count / stats.junction_count
    return stats


def extract_json_block(text: str, title: str) -> dict:
    marker = f"{title}:"
    idx = text.find(marker)
    if idx == -1:
        raise ValueError(f"Could not find {title!r} in command output")

    start = text.find("{", idx)
    if start == -1:
        raise ValueError(f"Could not find JSON start after {title!r}")

    depth = 0
    end = None
    for pos in range(start, len(text)):
        ch = text[pos]
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                end = pos + 1
                break

    if end is None:
        raise ValueError(f"Could not find JSON end for {title!r}")

    return json.loads(text[start:end])


def run_make_target(repo_root: Path, target: str) -> tuple[int, str]:
    proc = subprocess.run(
        ["make", target],
        cwd=repo_root,
        capture_output=True,
        text=True,
    )
    output = proc.stdout + proc.stderr
    return proc.returncode, output


def compare_stats(expected: IndexStats, actual: dict, rel_tol: float) -> list[str]:
    mismatches = []
    expected_dict = asdict(expected)
    for key, expected_value in expected_dict.items():
        if key not in actual:
            mismatches.append(f"missing key in actual summary: {key}")
            continue

        actual_value = actual[key]
        if isinstance(expected_value, float):
            if abs(expected_value - float(actual_value)) > rel_tol:
                mismatches.append(
                    f"{key}: expected {expected_value:.12f}, actual {float(actual_value):.12f}"
                )
        else:
            if expected_value != actual_value:
                mismatches.append(f"{key}: expected {expected_value}, actual {actual_value}")

    for key in actual:
        if key not in expected_dict:
            mismatches.append(f"unexpected key in actual summary: {key}")

    return mismatches


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Verify isomatch index summary against stats recomputed from the source GTF."
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Repository root containing the Makefile",
    )
    parser.add_argument(
        "--gtf",
        type=Path,
        default=Path("test/gencode.v49.basic.annotation.sorted.gtf.gz"),
        help="GTF used by the index command",
    )
    parser.add_argument(
        "--reffa",
        type=Path,
        default=Path("test/GRCh38.p14.allChr.fa"),
        help="Reference FASTA used for canonical junction checks",
    )
    parser.add_argument(
        "--make-target",
        default="index2",
        help="Make target to run for collecting the Rust summary",
    )
    parser.add_argument(
        "--float-tol",
        type=float,
        default=1e-12,
        help="Absolute tolerance for floating-point comparisons",
    )
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    gtf_path = (repo_root / args.gtf).resolve() if not args.gtf.is_absolute() else args.gtf
    reffa_path = (repo_root / args.reffa).resolve() if not args.reffa.is_absolute() else args.reffa

    expected = compute_expected_stats(gtf_path, reffa_path)
    expected_json = asdict(expected)

    print("Expected summary:")
    print(json.dumps(expected_json, indent=2, sort_keys=True))

    returncode, output = run_make_target(repo_root, args.make_target)
    print(f"\nCommand summary source: make {args.make_target}")
    print(f"Command exit code: {returncode}")

    if returncode != 0:
        sys.stderr.write(output)
        print("\nVerification failed because make returned a non-zero exit code.")
        return 1

    actual_json = extract_json_block(output, "Index summary")
    print("\nActual summary:")
    print(json.dumps(actual_json, indent=2, sort_keys=True))

    mismatches = compare_stats(expected, actual_json, args.float_tol)
    if mismatches:
        print("\nMismatch detected:")
        for mismatch in mismatches:
            print(f"- {mismatch}")
        return 1

    print("\nVerification passed: expected and actual index summaries match.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
