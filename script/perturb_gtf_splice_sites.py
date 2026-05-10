#!/usr/bin/env python3
"""
Randomly perturb splice sites for a fraction of transcripts in a GTF file.

The script performs two streaming passes over the input:

1. Collect eligible transcript IDs (multi-exon transcripts only).
2. Sample a fraction of those transcripts and write a new GTF with randomly
   shifted splice sites for the selected transcript blocks.

For each selected transcript, every splice junction is perturbed by shifting:
- the end coordinate of the upstream exon
- the start coordinate of the downstream exon

The transcript record receives a summary attribute, and the touched exon records
receive per-boundary shift attributes.
"""

from __future__ import annotations

import argparse
import gzip
import math
import random
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


TRANSCRIPT_ID_RE = re.compile(r'transcript_id\s+"([^"]+)"')


@dataclass
class GtfRecord:
    seqname: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    frame: str
    attrs: list[tuple[str, str]]


def open_text(path: Path, mode: str) -> TextIO:
    if "b" in mode:
        raise ValueError("binary mode is not supported")
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def extract_transcript_id(attr_text: str) -> str | None:
    match = TRANSCRIPT_ID_RE.search(attr_text)
    if match:
        return match.group(1)
    return None


def parse_attributes(attr_text: str) -> list[tuple[str, str]]:
    attrs: list[tuple[str, str]] = []
    for item in attr_text.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if " " not in item:
            attrs.append((item, ""))
            continue
        key, raw_value = item.split(" ", 1)
        value = raw_value.strip()
        if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
            value = value[1:-1]
        attrs.append((key, value))
    return attrs


def format_attributes(attrs: list[tuple[str, str]]) -> str:
    parts: list[str] = []
    for key, value in attrs:
        if value == "":
            parts.append(key)
        else:
            parts.append(f'{key} "{value}"')
    return "; ".join(parts) + ";"


def set_attribute(attrs: list[tuple[str, str]], key: str, value: str) -> None:
    for index, (existing_key, _) in enumerate(attrs):
        if existing_key == key:
            attrs[index] = (key, value)
            return
    attrs.append((key, value))


def parse_record(line: str) -> GtfRecord:
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 9:
        raise ValueError(f"expected 9 GTF columns, got {len(fields)}: {line.rstrip()}")
    return GtfRecord(
        seqname=fields[0],
        source=fields[1],
        feature=fields[2],
        start=int(fields[3]),
        end=int(fields[4]),
        score=fields[5],
        strand=fields[6],
        frame=fields[7],
        attrs=parse_attributes(fields[8]),
    )


def format_record(record: GtfRecord) -> str:
    return "\t".join(
        [
            record.seqname,
            record.source,
            record.feature,
            str(record.start),
            str(record.end),
            record.score,
            record.strand,
            record.frame,
            format_attributes(record.attrs),
        ]
    )


def iter_blocks(path: Path) -> Iterator[tuple[str, str | None, list[str]]]:
    pending_comments: list[str] = []
    current_tx_id: str | None = None
    current_block: list[str] = []

    with open_text(path, "rt") as handle:
        for raw_line in handle:
            if raw_line.startswith("#"):
                if current_block:
                    current_block.append(raw_line)
                else:
                    pending_comments.append(raw_line)
                continue

            fields = raw_line.rstrip("\n").split("\t", 8)
            if len(fields) != 9:
                raise ValueError(f"expected 9 GTF columns, got {len(fields)}: {raw_line.rstrip()}")

            feature = fields[2]
            tx_id = extract_transcript_id(fields[8])

            if feature == "transcript":
                if pending_comments:
                    yield ("comment", None, pending_comments)
                    pending_comments = []
                if current_block:
                    yield ("transcript", current_tx_id, current_block)
                current_tx_id = tx_id
                current_block = [raw_line]
                continue

            if current_block:
                current_block.append(raw_line)
            else:
                if pending_comments:
                    yield ("comment", None, pending_comments)
                    pending_comments = []
                yield ("other", tx_id, [raw_line])

        if pending_comments:
            yield ("comment", None, pending_comments)
        if current_block:
            yield ("transcript", current_tx_id, current_block)


def collect_eligible_transcripts(path: Path) -> list[str]:
    eligible_ids: list[str] = []
    for kind, tx_id, block_lines in iter_blocks(path):
        if kind != "transcript" or tx_id is None:
            continue
        exon_count = 0
        for raw_line in block_lines:
            if raw_line.startswith("#"):
                continue
            fields = raw_line.split("\t", 3)
            if len(fields) >= 3 and fields[2] == "exon":
                exon_count += 1
                if exon_count >= 2:
                    eligible_ids.append(tx_id)
                    break
    return eligible_ids


def make_output_path(base_output: Path, replicate_index: int, replicates: int) -> Path:
    if replicates == 1:
        return base_output

    name = base_output.name
    if name.endswith(".gtf.gz"):
        stem = name[: -len(".gtf.gz")]
        new_name = f"{stem}.rep{replicate_index}.gtf.gz"
    elif name.endswith(".gtf"):
        stem = name[: -len(".gtf")]
        new_name = f"{stem}.rep{replicate_index}.gtf"
    else:
        new_name = f"{name}.rep{replicate_index}.gtf"
    return base_output.with_name(new_name)


def choose_shift_pair(
    left_start: int,
    left_end: int,
    right_start: int,
    right_end: int,
    max_shift_bp: int,
    rng: random.Random,
) -> tuple[int, int]:
    end_delta_min = max(-max_shift_bp, left_start - left_end)
    end_delta_max = max_shift_bp
    start_delta_min = -max_shift_bp
    start_delta_max = min(max_shift_bp, right_end - right_start)
    gap = right_start - left_end
    valid_end_max = min(end_delta_max, start_delta_max - 2 + gap)

    if end_delta_min > valid_end_max:
        return 0, 0

    for require_both_nonzero in (True, False):
        for _ in range(64):
            end_delta = random_int_prefer_nonzero(end_delta_min, valid_end_max, rng)
            valid_start_min = max(start_delta_min, 2 - gap + end_delta)
            if valid_start_min > start_delta_max:
                continue
            start_delta = random_int_prefer_nonzero(valid_start_min, start_delta_max, rng)
            if require_both_nonzero and end_delta != 0 and start_delta != 0:
                return end_delta, start_delta
            if not require_both_nonzero and (end_delta != 0 or start_delta != 0):
                return end_delta, start_delta

    for end_delta in preferred_ints(end_delta_min, valid_end_max):
        valid_start_min = max(start_delta_min, 2 - gap + end_delta)
        if valid_start_min > start_delta_max:
            continue
        for start_delta in preferred_ints(valid_start_min, start_delta_max):
            if end_delta != 0 and start_delta != 0:
                return end_delta, start_delta
    for end_delta in preferred_ints(end_delta_min, valid_end_max):
        valid_start_min = max(start_delta_min, 2 - gap + end_delta)
        if valid_start_min > start_delta_max:
            continue
        for start_delta in preferred_ints(valid_start_min, start_delta_max):
            if end_delta != 0 or start_delta != 0:
                return end_delta, start_delta

    return 0, 0


def preferred_ints(low: int, high: int) -> list[int]:
    values: list[int] = []
    for candidate in (-1, 1, low, high, 0):
        if low <= candidate <= high and candidate not in values:
            values.append(candidate)
    return values


def random_int_prefer_nonzero(low: int, high: int, rng: random.Random) -> int:
    if low > high:
        raise ValueError(f"invalid integer range: {low}..{high}")
    if low == high:
        return low
    if 0 < low or 0 > high:
        return rng.randint(low, high)

    nonzero_count = (high - low + 1) - 1
    if nonzero_count <= 0:
        return 0

    offset = rng.randint(0, nonzero_count - 1)
    value = low + offset
    if value >= 0:
        value += 1
    return value


def perturb_transcript_block(block_lines: list[str], max_shift_bp: int, rng: random.Random) -> list[str]:
    records = [parse_record(line) for line in block_lines if not line.startswith("#")]
    comments = [line for line in block_lines if line.startswith("#")]

    transcript_record = next((record for record in records if record.feature == "transcript"), None)
    exon_records = [record for record in records if record.feature == "exon"]

    if transcript_record is None or len(exon_records) < 2:
        return block_lines

    original_starts = [record.start for record in exon_records]
    original_ends = [record.end for record in exon_records]
    start_shifts = [0] * len(exon_records)
    end_shifts = [0] * len(exon_records)
    junction_summaries: list[str] = []

    for junction_index in range(len(exon_records) - 1):
        left = exon_records[junction_index]
        right = exon_records[junction_index + 1]
        end_delta, start_delta = choose_shift_pair(
            left_start=left.start,
            left_end=left.end,
            right_start=right.start,
            right_end=right.end,
            max_shift_bp=max_shift_bp,
            rng=rng,
        )
        old_left_end = left.end
        old_right_start = right.start
        left.end += end_delta
        right.start += start_delta
        end_shifts[junction_index] += end_delta
        start_shifts[junction_index + 1] += start_delta
        junction_summaries.append(
            (
                f"J{junction_index + 1}:"
                f"left={old_left_end}->{left.end}({end_delta:+d}),"
                f"right={old_right_start}->{right.start}({start_delta:+d})"
            )
        )

    transcript_record.start = exon_records[0].start
    transcript_record.end = exon_records[-1].end
    set_attribute(transcript_record.attrs, "splice_site_perturbed", "true")
    set_attribute(transcript_record.attrs, "splice_junction_count", str(len(exon_records) - 1))
    set_attribute(transcript_record.attrs, "splice_junction_shifts", "|".join(junction_summaries))

    for exon_index, exon_record in enumerate(exon_records):
        if start_shifts[exon_index] != 0:
            set_attribute(exon_record.attrs, "splice_site_start_shift_bp", str(start_shifts[exon_index]))
        if end_shifts[exon_index] != 0:
            set_attribute(exon_record.attrs, "splice_site_end_shift_bp", str(end_shifts[exon_index]))
        if start_shifts[exon_index] != 0 or end_shifts[exon_index] != 0:
            set_attribute(exon_record.attrs, "original_exon_start", str(original_starts[exon_index]))
            set_attribute(exon_record.attrs, "original_exon_end", str(original_ends[exon_index]))

    comment_iter = iter(comments)
    formatted_records = [format_record(record) + "\n" for record in records]
    if not comments:
        return formatted_records

    output_lines: list[str] = []
    inserted_comments = False
    for line in formatted_records:
        output_lines.append(line)
        if not inserted_comments:
            output_lines.extend(comment_iter)
            inserted_comments = True
    output_lines.extend(comment_iter)
    return output_lines


def write_perturbed_gtf(
    input_path: Path,
    output_path: Path,
    selected_ids: set[str],
    max_shift_bp: int,
    rng: random.Random,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(output_path, "wt") as out_handle:
        for kind, tx_id, block_lines in iter_blocks(input_path):
            if kind == "comment" or kind == "other" or tx_id not in selected_ids:
                out_handle.writelines(block_lines)
                continue
            out_handle.writelines(perturb_transcript_block(block_lines, max_shift_bp=max_shift_bp, rng=rng))


def positive_fraction(value: str) -> float:
    parsed = float(value)
    if not 0.0 <= parsed <= 1.0:
        raise argparse.ArgumentTypeError("fraction must be between 0 and 1")
    return parsed


def non_negative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be >= 0")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Randomly perturb splice sites for a fraction of multi-exon "
            "transcripts in a GTF file."
        )
    )
    parser.add_argument("input_gtf", type=Path, help="Input GTF or GTF.GZ file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output GTF path or output prefix")
    parser.add_argument(
        "--fraction",
        type=positive_fraction,
        required=True,
        help="Fraction of eligible transcripts to perturb, between 0 and 1",
    )
    parser.add_argument(
        "--max-shift-bp",
        type=non_negative_int,
        required=True,
        help="Maximum absolute splice-site shift in base pairs",
    )
    parser.add_argument(
        "--replicates",
        type=non_negative_int,
        default=1,
        help="Number of independently perturbed output files to generate",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Random seed for reproducible transcript sampling and coordinate shifts",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.replicates < 1:
        raise SystemExit("--replicates must be at least 1")

    eligible_ids = collect_eligible_transcripts(args.input_gtf)
    eligible_count = len(eligible_ids)
    selected_count = min(eligible_count, int(math.floor(eligible_count * args.fraction + 0.5)))

    print(f"Eligible multi-exon transcripts: {eligible_count}")
    print(f"Selected transcripts per replicate: {selected_count}")

    for replicate_index in range(1, args.replicates + 1):
        replicate_seed = args.seed + replicate_index - 1
        rng = random.Random(replicate_seed)
        selected_ids = set(rng.sample(eligible_ids, selected_count)) if selected_count > 0 else set()
        output_path = make_output_path(args.output, replicate_index, args.replicates)
        write_perturbed_gtf(
            input_path=args.input_gtf,
            output_path=output_path,
            selected_ids=selected_ids,
            max_shift_bp=args.max_shift_bp,
            rng=rng,
        )
        print(f"[replicate {replicate_index}] seed={replicate_seed} wrote {output_path}")


if __name__ == "__main__":
    main()
