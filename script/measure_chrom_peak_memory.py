#!/usr/bin/env python3
"""
CODEX:
Estimate peak memory while reading one chromosome with the current reader model.

This script does not use `cargo`, so it stays usable even when the Rust crate
has unrelated compile errors. It parses the `.isomx` header/directory directly
and estimates memory based on how `IndexReader` currently loads sections:

- `TxBase` records are streamed one by one from disk.
- `JunctionPool::load_range()` reads a raw byte buffer, then materializes a
  decoded `Vec<u32>`, so peak during that call is roughly `2 * section_len`.
- `StringPool::load_range()` keeps the raw bytes as-is, so peak is roughly
  `section_len`.
- `SpliceSitePool::load_range()` reads a raw byte buffer, then materializes a
  decoded `Vec<SpliceSitePair>`, so peak during that call is roughly
  `2 * section_len`.
"""

from __future__ import annotations

import argparse
import struct
from dataclasses import dataclass
from pathlib import Path


HEADER_SIZE = 4096
DIR_ENTRY_SIZE = 42
TX_BASE_DISK_SIZE = 76


@dataclass
class ChromEntry:
    chrom_id: int
    chrom_name: str
    tx_count: int
    tx_offset: int
    junction_offset: int
    junction_len: int
    string_offset: int
    string_len: int
    splice_offset: int
    splice_len: int


def fmt_bytes(n: int) -> str:
    mib = n / 1024 / 1024
    return f"{n:,} bytes ({mib:.2f} MiB)"


def parse_index(path: Path) -> list[ChromEntry]:
    with path.open("rb") as fh:
        header = fh.read(HEADER_SIZE)
        chrom_count = struct.unpack("<I", header[13:17])[0]
        chrom_name_table_len = struct.unpack("<I", header[41:45])[0]

        raw_entries = [struct.unpack("<HIIIIIIIIII", fh.read(DIR_ENTRY_SIZE)) for _ in range(chrom_count)]
        chrom_name_table = fh.read(chrom_name_table_len)

    entries: list[ChromEntry] = []
    for raw in raw_entries:
        (
            chrom_id,
            chrom_name_offset,
            chrom_name_len,
            tx_count,
            tx_offset,
            junction_offset,
            junction_len,
            string_offset,
            string_len,
            splice_offset,
            splice_len,
        ) = raw
        chrom_name = chrom_name_table[chrom_name_offset : chrom_name_offset + chrom_name_len].decode()
        entries.append(
            ChromEntry(
                chrom_id=chrom_id,
                chrom_name=chrom_name,
                tx_count=tx_count,
                tx_offset=tx_offset,
                junction_offset=junction_offset,
                junction_len=junction_len,
                string_offset=string_offset,
                string_len=string_len,
                splice_offset=splice_offset,
                splice_len=splice_len,
            )
        )
    return entries


def find_chrom(entries: list[ChromEntry], chrom: str) -> ChromEntry:
    normalized = chrom.lower()
    for entry in entries:
        if entry.chrom_name.lower() == normalized:
            return entry
    raise SystemExit(f"chromosome not found in index: {chrom}")


def print_report(entry: ChromEntry) -> None:
    tx_streaming_working_set = TX_BASE_DISK_SIZE

    # Current implementation behavior:
    # - load_junction_pool: buf(len) + Vec<u32>(len)
    # - load_string_pool:   Vec<u8>(len)
    # - load_splice_pool:   buf(len) + Vec<u8>(len)
    after_junction_resident = entry.junction_len
    junction_call_peak = entry.junction_len * 2

    after_string_resident = after_junction_resident + entry.string_len
    string_call_peak = after_junction_resident + entry.string_len

    after_splice_resident = after_string_resident + entry.splice_len
    splice_call_peak = after_string_resident + entry.splice_len * 2

    estimated_peak = max(
        tx_streaming_working_set,
        junction_call_peak,
        string_call_peak,
        splice_call_peak,
    )

    print(f"Chromosome: {entry.chrom_name} (id={entry.chrom_id})")
    print(f"Transcript count: {entry.tx_count:,}")
    print(f"Tx bytes on disk: {fmt_bytes(entry.tx_count * TX_BASE_DISK_SIZE)}")
    print(f"Junction pool bytes: {fmt_bytes(entry.junction_len)}")
    print(f"String pool bytes: {fmt_bytes(entry.string_len)}")
    print(f"Splice-site pool bytes: {fmt_bytes(entry.splice_len)}")
    print()
    print("Estimated resident memory by stage:")
    print(f"- streaming TxBase iteration only: {fmt_bytes(tx_streaming_working_set)}")
    print(f"- after load_junction_pool(): {fmt_bytes(after_junction_resident)}")
    print(f"- after load_string_pool(): {fmt_bytes(after_string_resident)}")
    print(f"- after load_splice_site_pool(): {fmt_bytes(after_splice_resident)}")
    print()
    print("Estimated peak during load calls:")
    print(f"- inside load_junction_pool(): {fmt_bytes(junction_call_peak)}")
    print(f"- inside load_string_pool(): {fmt_bytes(string_call_peak)}")
    print(f"- inside load_splice_site_pool(): {fmt_bytes(splice_call_peak)}")
    print()
    print(f"Estimated peak for current read path: {fmt_bytes(estimated_peak)}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Estimate peak memory for reading one chromosome from an .isomx index")
    parser.add_argument("index_path", type=Path, help="Path to the .isomx file")
    parser.add_argument("--chrom", default="chr1", help="Chromosome name to inspect, default: chr1")
    args = parser.parse_args()

    entries = parse_index(args.index_path)
    entry = find_chrom(entries, args.chrom)
    print_report(entry)


if __name__ == "__main__":
    main()
