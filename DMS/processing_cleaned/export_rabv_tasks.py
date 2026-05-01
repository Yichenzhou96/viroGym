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


def build_rabv_rows():
    reference_sequence = parse_fasta_sequence(raw_csv_path("RABV", "reference.fasta"))
    source_path = raw_csv_path("RABV", "all_antibodies_and_cell_entry.csv")

    antibody_rows = []
    cell_entry_rows = []
    with open(source_path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            site = int(row["sequential_site"])
            mutation = f"{row['wildtype']}{site}{row['mutant']}"
            mutated_sequence = mutate_sequence(reference_sequence, site, row["mutant"])

            antibody_score = row["monoclonal antibodies escape"]
            if antibody_score != "":
                antibody_rows.append(
                    {
                        "mutation": mutation,
                        "mutated_sequence": mutated_sequence,
                        "DMS_score": float(antibody_score),
                    }
                )

            cell_entry_score = row["cell entry"]
            if cell_entry_score != "":
                cell_entry_rows.append(
                    {
                        "mutation": mutation,
                        "mutated_sequence": mutated_sequence,
                        "DMS_score": float(cell_entry_score),
                    }
                )

    return {
        "RABV_antibody_escape.csv": antibody_rows,
        "RABV_cell_entry.csv": cell_entry_rows,
    }


def main():
    parser = argparse.ArgumentParser(description="Export canonical RABV task CSVs.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path.cwd(),
        help="Directory where canonical RABV CSVs will be written.",
    )
    args = parser.parse_args()

    outputs = build_rabv_rows()
    for filename, rows in outputs.items():
        destination = args.outdir / filename
        write_rows(destination, rows)
        print(f"Wrote {len(rows)} rows to {destination}")


if __name__ == "__main__":
    main()
