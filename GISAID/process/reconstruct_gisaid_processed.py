#!/usr/bin/env python3
"""Rebuild base GISAID processed CSVs from reference sequences and yearly counts.

This reconstructs the non-model columns for:

* GISAID_COV2_mapped_2021_Reference_processed.csv
* GISAID_COV2_mapped_2022_Alpha_processed.csv
* GISAID_COV2_mapped_2023_BA1_processed.csv
* GISAID_COV2_mapped_2024_XBB15_processed.csv
* GISAID_COV2_mapped_2025_KP3_processed.csv

The model-scoring columns are appended later by the in-place ESM and ProGen2
scripts in this directory.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


AMINO_ACIDS = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "Q",
    "E",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
]

WINDOWS = [
    {
        "output_name": "2021_Reference",
        "reference_name": "Reference",
        "count_t_json": "data_2020.json",
        "count_t1_json": "data_2021.json",
    },
    {
        "output_name": "2022_Alpha",
        "reference_name": "20I_Alpha_V1",
        "count_t_json": "data_2021.json",
        "count_t1_json": "data_2022.json",
    },
    {
        "output_name": "2023_BA1",
        "reference_name": "21K_Omicron_BA.1",
        "count_t_json": "data_2022.json",
        "count_t1_json": "data_2023.json",
    },
    {
        "output_name": "2024_XBB15",
        "reference_name": "XBB.1.5",
        "count_t_json": "data_2023.json",
        "count_t1_json": "data_2024.json",
    },
    {
        "output_name": "2025_KP3",
        "reference_name": "KP.3",
        "count_t_json": "data_2024.json",
        "count_t1_json": "data_2025.json",
    },
]

FIELDNAMES = [
    "i",
    "wt",
    "mut",
    "mutated_sequence",
    "mutation",
    "count_T",
    "count_T1",
    "freq_T",
    "freq_T1",
    "freq_delta",
    "fold_change",
    "is_novel",
    "is_emergent",
]


def load_reference_sequences(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {row["name"]: row["sequence"] for row in rows}


def load_counts(path: Path) -> tuple[dict[str, int], int]:
    with path.open() as handle:
        counts = json.load(handle)
    return counts, sum(int(value) for value in counts.values())


def count_lookup_key(wuhan_wt: str, position: int, target: str) -> str | None:
    if target == wuhan_wt:
        return None
    if target == "-":
        return f"{wuhan_wt}{position}del"
    return f"{wuhan_wt}{position}{target}"


def compute_fold_change(freq_t: float, freq_t1: float, count_t: int, count_t1: int) -> float:
    if count_t == 0:
        return math.inf if count_t1 > 0 else 0.0
    return round(freq_t1 / freq_t, 2)


def is_emergent(count_t: int, fold_change: float, freq_delta: float) -> int:
    if freq_delta >= 1e-4:
        return 1
    if count_t > 0 and fold_change >= 10.0:
        return 1
    return 0


def format_frequency(value: float) -> str:
    if value == 0.0:
        return "0.0"
    return format(value, ".7g")


def iter_rows(
    reference_sequence: str,
    wuhan_sequence: str,
    counts_t: dict[str, int],
    counts_t1: dict[str, int],
    total_t: int,
    total_t1: int,
):
    if len(reference_sequence) != len(wuhan_sequence):
        raise ValueError("Reference and Wuhan sequences must have the same length")

    for position, (reference_char, wuhan_char) in enumerate(zip(reference_sequence, wuhan_sequence), start=1):
        targets = [aa for aa in AMINO_ACIDS if aa != reference_char]
        if reference_char != "-":
            targets.append("-")

        for target in targets:
            mutation_label = f"{reference_char}{position}{target}"
            key = count_lookup_key(wuhan_char, position, target)
            count_t = int(counts_t.get(key, 0)) if key else 0
            count_t1 = int(counts_t1.get(key, 0)) if key else 0

            freq_t = count_t / total_t if total_t else 0.0
            freq_t1 = count_t1 / total_t1 if total_t1 else 0.0
            freq_delta = freq_t1 - freq_t
            fold_change = compute_fold_change(freq_t, freq_t1, count_t, count_t1)
            novel = int(count_t == 0 and count_t1 > 0)
            emergent = is_emergent(count_t, fold_change, freq_delta)

            seq_chars = list(reference_sequence)
            seq_chars[position - 1] = target

            yield {
                "i": position,
                "wt": wuhan_char,
                "mut": target,
                "mutated_sequence": "".join(seq_chars),
                "mutation": mutation_label,
                "count_T": count_t,
                "count_T1": count_t1,
                "freq_T": format_frequency(freq_t),
                "freq_T1": format_frequency(freq_t1),
                "freq_delta": format_frequency(freq_delta),
                "fold_change": fold_change,
                "is_novel": novel,
                "is_emergent": emergent,
            }


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def build_outputs(gisaid_dir: Path, workspace_root: Path, output_dir: Path) -> None:
    reference_sequences = load_reference_sequences(gisaid_dir / "reference_sequence.csv")
    wuhan_sequence = reference_sequences["Reference"]

    for window in WINDOWS:
        counts_t, total_t = load_counts(workspace_root / window["count_t_json"])
        counts_t1, total_t1 = load_counts(workspace_root / window["count_t1_json"])
        rows = list(
            iter_rows(
                reference_sequence=reference_sequences[window["reference_name"]],
                wuhan_sequence=wuhan_sequence,
                counts_t=counts_t,
                counts_t1=counts_t1,
                total_t=total_t,
                total_t1=total_t1,
            )
        )
        output_path = output_dir / f"GISAID_COV2_mapped_{window['output_name']}_processed.csv"
        write_rows(output_path, rows)
        print(f"wrote {output_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gisaid-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory containing reference_sequence.csv.",
    )
    parser.add_argument(
        "--workspace-root",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
        help="Directory containing data_*.json source files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory to write reconstructed processed CSVs into.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    build_outputs(
        gisaid_dir=args.gisaid_dir.resolve(),
        workspace_root=args.workspace_root.resolve(),
        output_dir=args.output_dir.resolve(),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
