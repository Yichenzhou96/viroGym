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
            reference_site = row["reference_site"].strip()
            site_map[reference_site] = {
                "sequential_site": int(float(row["sequential_site"])),
                "sequential_wt": row["sequential_wt"].strip(),
            }
    return site_map


def group_rows(rows):
    grouped = {}
    for row in rows:
        key = (row["mutation"], row["mutated_sequence"])
        grouped.setdefault(key, []).append(float(row["DMS_score"]))

    output_rows = []
    for (mutation, mutated_sequence), scores in grouped.items():
        output_rows.append(
            {
                "mutation": mutation,
                "mutated_sequence": mutated_sequence,
                "DMS_score": mean(scores),
            }
        )
    return output_rows


def merge_mean(left_rows, right_rows):
    right_by_key = {
        (row["mutation"], row["mutated_sequence"]): row for row in right_rows
    }
    merged_rows = []
    seen = set()

    for row in left_rows:
        key = (row["mutation"], row["mutated_sequence"])
        seen.add(key)
        if key in right_by_key:
            score = mean([row["DMS_score"], right_by_key[key]["DMS_score"]])
        else:
            score = row["DMS_score"]
        merged_rows.append(
            {
                "mutation": row["mutation"],
                "mutated_sequence": row["mutated_sequence"],
                "DMS_score": score,
            }
        )

    for row in right_rows:
        key = (row["mutation"], row["mutated_sequence"])
        if key in seen:
            continue
        merged_rows.append(
            {
                "mutation": row["mutation"],
                "mutated_sequence": row["mutated_sequence"],
                "DMS_score": row["DMS_score"],
            }
        )

    return merged_rows


def load_antibody_rows(path, site_map, reference_sequence):
    raw_rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            site_info = site_map[row["site"].strip()]
            if site_info["sequential_wt"] == "*":
                continue
            site = site_info["sequential_site"]
            mutant = row["mutant"].strip()
            mutation = f"{site_info['sequential_wt']}{site}{mutant}"
            raw_rows.append(
                {
                    "mutation": mutation,
                    "mutated_sequence": mutate_sequence(reference_sequence, site, mutant),
                    "DMS_score": float(row["escape_mean"]),
                }
            )
    return group_rows(raw_rows)


def load_cell_entry_rows(path, site_map, reference_sequence):
    raw_rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            site_info = site_map[row["site"].strip()]
            site = site_info["sequential_site"]
            sequential_wt = site_info["sequential_wt"]
            mutant = row["mutant"].strip()
            if mutant == "*" or sequential_wt == "*" or mutant == sequential_wt:
                continue
            mutation = f"{sequential_wt}{site}{mutant}"
            raw_rows.append(
                {
                    "mutation": mutation,
                    "mutated_sequence": mutate_sequence(reference_sequence, site, mutant),
                    "DMS_score": float(row["effect"]),
                }
            )
    return group_rows(raw_rows)


def drop_stop_rows(rows):
    output_rows = [row for row in rows if not row["mutation"].endswith("*")]
    output_rows.sort(key=lambda row: (row["mutation"], row["mutated_sequence"]))
    return output_rows


def build_bf520_antibody_rows():
    reference_sequence = parse_fasta_sequence(
        raw_csv_path("HIV", "mAb", "BF520", "reference.fasta")
    ).rstrip("*")
    site_map = load_site_map(raw_csv_path("HIV", "mAb", "BF520", "site_numbering_map.csv"))

    merged_rows = merge_mean(
        load_antibody_rows(raw_csv_path("HIV", "mAb", "BF520", "BF520_Env_3BNC117_mut_effects_filtered.csv"), site_map, reference_sequence),
        load_antibody_rows(raw_csv_path("HIV", "mAb", "BF520", "3BNC117_avg.csv"), site_map, reference_sequence),
    )
    merged_rows = merge_mean(
        merged_rows,
        load_antibody_rows(raw_csv_path("HIV", "mAb", "BF520", "PGT151_avg.csv"), site_map, reference_sequence),
    )
    merged_rows = merge_mean(
        merged_rows,
        load_antibody_rows(raw_csv_path("HIV", "mAb", "BF520", "1-18_avg.csv"), site_map, reference_sequence),
    )
    merged_rows = merge_mean(
        merged_rows,
        load_antibody_rows(raw_csv_path("HIV", "mAb", "BF520", "BF520_Env_10-1074_mut_effects_filtered.csv"), site_map, reference_sequence),
    )
    return drop_stop_rows(merged_rows)


def build_tro11_antibody_rows():
    reference_sequence = parse_fasta_sequence(
        raw_csv_path("HIV", "mAb", "TRO11", "reference.fasta")
    ).rstrip("*")
    site_map = load_site_map(raw_csv_path("HIV", "mAb", "TRO11", "site_numbering_map.csv"))
    merged_rows = merge_mean(
        load_antibody_rows(raw_csv_path("HIV", "mAb", "TRO11", "TRO11_Env_3BNC117_mut_effects_filtered.csv"), site_map, reference_sequence),
        load_antibody_rows(raw_csv_path("HIV", "mAb", "TRO11", "TRO11_Env_10-1074_mut_effects_filtered.csv"), site_map, reference_sequence),
    )
    return drop_stop_rows(merged_rows)


def build_hiv_mab_rows():
    bf520_reference_sequence = parse_fasta_sequence(
        raw_csv_path("HIV", "mAb", "BF520", "reference.fasta")
    ).rstrip("*")
    bf520_site_map = load_site_map(raw_csv_path("HIV", "mAb", "BF520", "site_numbering_map.csv"))

    tro11_reference_sequence = parse_fasta_sequence(
        raw_csv_path("HIV", "mAb", "TRO11", "reference.fasta")
    ).rstrip("*")
    tro11_site_map = load_site_map(raw_csv_path("HIV", "mAb", "TRO11", "site_numbering_map.csv"))

    return {
        "HIV_antibody_escape_B520.csv": build_bf520_antibody_rows(),
        "HIV_antibody_escape_TRO11.csv": build_tro11_antibody_rows(),
        "HIV_cell_entry_B520.csv": load_cell_entry_rows(
            raw_csv_path("HIV", "mAb", "BF520", "BF520_Env_TZM-bl_entry_func_effects_filtered.csv"),
            bf520_site_map,
            bf520_reference_sequence,
        ),
        "HIV_cell_entry_TRO11.csv": load_cell_entry_rows(
            raw_csv_path("HIV", "mAb", "TRO11", "TRO11_Env_TZM-bl_entry_func_effects_filtered.csv"),
            tro11_site_map,
            tro11_reference_sequence,
        ),
    }


def main():
    parser = argparse.ArgumentParser(description="Export canonical HIV mAb task CSVs.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path.cwd(),
        help="Directory where canonical HIV mAb CSVs will be written.",
    )
    args = parser.parse_args()

    outputs = build_hiv_mab_rows()
    for filename, rows in outputs.items():
        destination = args.outdir / filename
        write_rows(destination, rows)
        print(f"Wrote {len(rows)} rows to {destination}")


if __name__ == "__main__":
    main()
