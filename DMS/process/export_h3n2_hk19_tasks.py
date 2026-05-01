#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
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
            reference_site = row["reference_site"]
            if reference_site == "":
                continue
            site_map[int(reference_site)] = (int(row["sequential_site"]), row["sequential_wt"])
    return site_map


def build_hk19_task_rows(score_column):
    reference_sequence = parse_fasta_sequence(
        raw_csv_path("influenza", "full_hk19", "hk19_reference.fasta")
    )
    site_map = load_site_map(raw_csv_path("influenza", "full_hk19", "site_map.csv"))

    scores_by_mutation = defaultdict(list)
    sequences_by_mutation = {}
    with open(raw_csv_path("influenza", "full_hk19", "full_hk19_escape_scores.csv"), newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            mutant = row["mutant"]
            if mutant in {"", "-"} or row[score_column] == "":
                continue
            raw_site = int(row["site"])
            if raw_site not in site_map:
                continue

            sequential_site, sequential_wildtype = site_map[raw_site]
            mutation = f"{sequential_wildtype}{sequential_site}{mutant}"
            sequences_by_mutation[mutation] = mutate_sequence(
                reference_sequence, sequential_site, mutant
            )
            scores_by_mutation[mutation].append(row[score_column])

    rows = []
    for mutation in sorted(scores_by_mutation):
        rows.append(
            {
                "mutation": mutation,
                "mutated_sequence": sequences_by_mutation[mutation],
                "DMS_score": mean(scores_by_mutation[mutation]),
            }
        )
    return rows


def build_hk19_rows():
    return {
        "FLU_sera_escape_H3N2_HK19.csv": build_hk19_task_rows("escape_mean"),
        "FLU_functional_effect_H3N2_HK19.csv": build_hk19_task_rows("functional effect"),
    }


def main():
    parser = argparse.ArgumentParser(description="Export canonical H3N2 HK19 influenza task CSVs.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path.cwd(),
        help="Directory where canonical H3N2 HK19 CSVs will be written.",
    )
    args = parser.parse_args()

    outputs = build_hk19_rows()
    for filename, rows in outputs.items():
        destination = args.outdir / filename
        write_rows(destination, rows)
        print(f"Wrote {len(rows)} rows to {destination}")


if __name__ == "__main__":
    main()
