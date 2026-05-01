#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

from processing_cleaned.common import (
    mutate_sequence,
    parse_fasta_sequence,
    raw_csv_path,
    write_rows,
)


COLUMN_TO_OUTPUT = {
    "mean antibodies": "SARS_antibody_escape_KP311.csv",
    "ACE2 binding": "SARS_binding_KP311.csv",
    "spike mediated entry": "SARS_cell_entry_KP311.csv",
}


def build_kp311_rows():
    outputs = {filename: [] for filename in COLUMN_TO_OUTPUT.values()}
    source_path = raw_csv_path("SARS-CoV-2", "KP311", "antibody_escape.csv")
    reference_sequence = parse_fasta_sequence(
        raw_csv_path("SARS-CoV-2", "KP311", "reference.fasta")
    ).rstrip("*")
    with open(source_path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            mutant = row["mutant"]
            if mutant in {"", "-"}:
                continue
            site = int(row["sequential_site"])
            mutation = f"{row['wildtype']}{site}{mutant}"
            mutated_sequence = mutate_sequence(reference_sequence, site, mutant)
            for column, filename in COLUMN_TO_OUTPUT.items():
                score = row[column]
                if score == "":
                    continue
                outputs[filename].append(
                    {
                        "mutation": mutation,
                        "mutated_sequence": mutated_sequence,
                        "DMS_score": float(score),
                    }
                )
    return outputs


def main():
    parser = argparse.ArgumentParser(description="Export canonical KP.3.1.1 SARS-CoV-2 task CSVs.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path.cwd(),
        help="Directory where canonical KP.3.1.1 CSVs will be written.",
    )
    args = parser.parse_args()

    outputs = build_kp311_rows()
    for filename, rows in outputs.items():
        destination = args.outdir / filename
        write_rows(destination, rows)
        print(f"Wrote {len(rows)} rows to {destination}")


if __name__ == "__main__":
    main()
