#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

from processing_cleaned.common import (
    mean,
    mutate_sequence,
    parse_fasta_sequence,
    raw_csv_path,
    write_rows,
)


def load_site_map(path):
    site_map = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            site_map[row["reference_site"].strip()] = int(float(row["sequential_site"]))
    return site_map


def load_serum_rows(path, site_map, ref_sequence):
    values = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            reference_site = row["site"].strip()
            site = site_map[reference_site]
            mutation = f"{row['wildtype']}{site}{row['mutant']}"
            mutated_sequence = mutate_sequence(ref_sequence, site, row["mutant"])
            key = (mutation, mutated_sequence)
            values.setdefault(key, []).append(float(row["escape_mean"]))
    return {key: mean(scores) for key, scores in values.items()}


def build_hiv_b520_rows():
    ref_sequence = parse_fasta_sequence(raw_csv_path("HIV", "pAb", "BF520", "reference.fasta"))
    site_map = load_site_map(raw_csv_path("HIV", "pAb", "BF520", "site_numbering_map.csv"))

    serum_513 = load_serum_rows(raw_csv_path("HIV", "pAb", "BF520", "IDC513_avg.csv"), site_map, ref_sequence)
    serum_561 = load_serum_rows(raw_csv_path("HIV", "pAb", "BF520", "IDC561_avg.csv"), site_map, ref_sequence)
    serum_033 = load_serum_rows(raw_csv_path("HIV", "pAb", "BF520", "IDF033_avg.csv"), site_map, ref_sequence)

    output_rows = []
    for mutation, mutated_sequence in serum_513:
        if (mutation, mutated_sequence) not in serum_561:
            continue
        if (mutation, mutated_sequence) not in serum_033:
            continue
        score = mean(
            [
                mean(
                    [
                        serum_513[(mutation, mutated_sequence)],
                        serum_561[(mutation, mutated_sequence)],
                    ]
                ),
                serum_033[(mutation, mutated_sequence)],
            ]
        )
        output_rows.append(
            {
                "mutation": mutation,
                "mutated_sequence": mutated_sequence,
                "DMS_score": score,
            }
        )
    return output_rows


def main():
    parser = argparse.ArgumentParser(description="Export canonical HIV_sera_escape_B520.csv rows.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path.cwd() / "HIV_sera_escape_B520.csv",
        help="Destination CSV path.",
    )
    args = parser.parse_args()

    rows = build_hiv_b520_rows()
    write_rows(args.output, rows)
    print(f"Wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
