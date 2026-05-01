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


def build_nipah_rows():
    reference_sequence = parse_fasta_sequence(raw_csv_path("Nipah", "reference.fasta")).rstrip("*")
    source_path = raw_csv_path("Nipah", "summary.csv")

    antibody_rows = []
    binding_rows = []
    cell_entry_rows = []
    with open(source_path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            site = int(row["sequential_site"])
            mutation = f"{row['wildtype']}{site}{row['mutant']}"
            mutated_sequence = mutate_sequence(reference_sequence, site, row["mutant"])

            antibody_score = row["antibodies escape"]
            if antibody_score != "":
                antibody_rows.append(
                    {
                        "mutation": mutation,
                        "mutated_sequence": mutated_sequence,
                        "DMS_score": float(antibody_score),
                    }
                )

            binding_values = [row["EphrinB2 binding"], row["EphrinB3 binding"]]
            if any(value != "" for value in binding_values):
                binding_rows.append(
                    {
                        "mutation": mutation,
                        "mutated_sequence": mutated_sequence,
                        "DMS_score": mean(binding_values),
                    }
                )

            entry_values = [row["CHO_bEFNB2_entry"], row["CHO_bEFNB3_entry"]]
            if any(value != "" for value in entry_values):
                cell_entry_rows.append(
                    {
                        "mutation": mutation,
                        "mutated_sequence": mutated_sequence,
                        "DMS_score": mean(entry_values),
                    }
                )

    return {
        "NIPAH_antibody_escape.csv": antibody_rows,
        "NIPAH_binding.csv": binding_rows,
        "NIPAH_cell_entry.csv": cell_entry_rows,
    }


def main():
    parser = argparse.ArgumentParser(description="Export canonical Nipah task CSVs.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path.cwd(),
        help="Directory where canonical Nipah CSVs will be written.",
    )
    args = parser.parse_args()

    outputs = build_nipah_rows()
    for filename, rows in outputs.items():
        destination = args.outdir / filename
        write_rows(destination, rows)
        print(f"Wrote {len(rows)} rows to {destination}")


if __name__ == "__main__":
    main()
