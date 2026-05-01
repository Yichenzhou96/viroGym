#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
from pathlib import Path

from processing_cleaned.common import (
    exact_mean,
    mutate_sequence,
    parse_fasta_sequence,
    raw_csv_path,
    write_rows,
)


ANTIBODY_SOURCES = [
    ("121F_filtered_mut_effect.csv", "121F"),
    ("256A_filtered_mut_effect.csv", "256A"),
    ("2510C_filtered_mut_effect.csv", "2510C"),
    ("372D_filtered_mut_effect.csv", "372D"),
    ("377H_filtered_mut_effect.csv", "377H"),
    ("89F_filtered_mut_effect.csv", "89F"),
]


def _iter_filtered_antibody_rows(reference_sequence):
    for filename, label in ANTIBODY_SOURCES:
        with open(raw_csv_path("Lassa", filename), newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                if row["poor_cell_entry"] != "False" or row["escape"] == "":
                    continue
                site = int(row["site"])
                mutation = f"{row['wildtype']}{site}{row['mutant']}"
                yield label, {
                    "mutation": mutation,
                    "mutated_sequence": mutate_sequence(reference_sequence, site, row["mutant"]),
                    "DMS_score": row["escape"],
                }


def _legacy_aliases(label):
    return [
        f"LASV_{label}_escape.csv",
        f"LASV_antibody_escape_{label}.csv",
    ]


def build_lassa_rows():
    reference_sequence = parse_fasta_sequence(raw_csv_path("Lassa", "reference.fasta")).rstrip("*")

    antibody_scores = defaultdict(list)
    for _, row in _iter_filtered_antibody_rows(reference_sequence):
        key = (row["mutation"], row["mutated_sequence"])
        antibody_scores[key].append(row["DMS_score"])

    antibody_rows = [
        {
            "mutation": mutation,
            "mutated_sequence": mutated_sequence,
            "DMS_score": exact_mean(scores),
        }
        for (mutation, mutated_sequence), scores in antibody_scores.items()
    ]
    antibody_rows.sort(key=lambda row: (row["mutation"], row["mutated_sequence"]))

    cell_entry_rows = []
    with open(raw_csv_path("Lassa", "293T_filtered_func_effects.csv"), newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            site = int(row["site"])
            wildtype = row["wildtype"]
            mutant = row["mutant"]
            if mutant in {"*", wildtype}:
                continue
            cell_entry_rows.append(
                {
                    "mutation": f"{wildtype}{site}{mutant}",
                    "mutated_sequence": mutate_sequence(reference_sequence, site, mutant),
                    "DMS_score": float(row["effect"]),
                }
            )

    return {
        "LASV_antibody_escape.csv": antibody_rows,
        "LASV_cell_entry.csv": cell_entry_rows,
    }


def build_legacy_lassa_rows():
    reference_sequence = parse_fasta_sequence(raw_csv_path("Lassa", "reference.fasta"))

    outputs = {}
    for label, row in _iter_filtered_antibody_rows(reference_sequence):
        for filename in _legacy_aliases(label):
            outputs.setdefault(filename, []).append(row)
    return outputs


def main():
    parser = argparse.ArgumentParser(description="Export canonical LASV task CSVs.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path.cwd(),
        help="Directory where canonical LASV CSVs will be written.",
    )
    parser.add_argument(
        "--include-legacy-files",
        action="store_true",
        help="Also write the legacy per-antibody LASV CSV aliases under the output directory.",
    )
    args = parser.parse_args()

    outputs = build_lassa_rows()
    if args.include_legacy_files:
        outputs.update(build_legacy_lassa_rows())
    for filename, rows in outputs.items():
        destination = args.outdir / filename
        write_rows(destination, rows)
        print(f"Wrote {len(rows)} rows to {destination}")


if __name__ == "__main__":
    main()
