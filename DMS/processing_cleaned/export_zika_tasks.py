#!/usr/bin/env python3

import argparse
import csv
from collections import OrderedDict
from pathlib import Path

from processing_cleaned.common import (
    infer_reference_sequence_from_task_rows,
    mean,
    mutate_sequence,
    raw_csv_path,
    write_rows,
)


ANTIBODY_SOURCES = [
    "antibody_EDE1-C10_median.csv",
    "antibody_EDE1-C8_median.csv",
    "antibody_MZ4_median.csv",
    "antibody_SIgN-3C_median.csv",
    "antibody_ZV-67_median.csv",
    "summary_ZKA64-medianmutdiffsel.csv",
    "summary_ZKA185-medianmutdiffsel.csv",
]


REFERENCE_TASK_PATH = raw_csv_path("ZIKA", "ZIKV_antibody_escapge_EDE1-C10.csv")


def build_zika_viral_growth_rows(reference_sequence):
    rows = []
    with open(raw_csv_path("ZIKA", "unscaled_muteffects.csv"), newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            site = int(row["site"])
            mutation = row["mutation"]
            mutated_sequence = mutate_sequence(reference_sequence, site, row["mutant"])
            rows.append(
                {
                    "mutation": mutation,
                    "mutated_sequence": mutated_sequence,
                    "DMS_score": float(row["effect"]),
                }
            )
    return rows


def load_antibody_table(path):
    values = OrderedDict()
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = (row["wildtype"], int(row["site"]), row["mutation"])
            values[key] = row["mutdiffsel"]
    return values


def load_antibody_tables():
    tables = [load_antibody_table(raw_csv_path("ZIKA", filename)) for filename in ANTIBODY_SOURCES]
    expected_keys = set(tables[0])
    for filename, table in zip(ANTIBODY_SOURCES[1:], tables[1:]):
        if set(table) != expected_keys:
            raise ValueError(f"ZIKV antibody table keys do not match baseline table: {filename}")
    return tables


def build_zika_antibody_rows(reference_sequence):
    tables = load_antibody_tables()
    rows = []
    for (wildtype, site, mutant), _ in tables[0].items():
        scores = [table[(wildtype, site, mutant)] for table in tables]
        if all(score == "" for score in scores):
            continue
        mutation = f"{wildtype}{site}{mutant}"
        mutated_sequence = mutate_sequence(reference_sequence, site, mutant)
        rows.append(
            {
                "mutation": mutation,
                "mutated_sequence": mutated_sequence,
                "DMS_score": mean(scores),
            }
        )
    rows.sort(key=lambda row: (row["mutation"], row["mutated_sequence"]))
    return rows


def build_zika_rows():
    reference_sequence = infer_reference_sequence_from_task_rows(REFERENCE_TASK_PATH)
    return {
        "ZIKV_antibody_escape.csv": build_zika_antibody_rows(reference_sequence),
        "ZIKV_viral_growth.csv": build_zika_viral_growth_rows(reference_sequence),
    }


def main():
    parser = argparse.ArgumentParser(description="Export canonical ZIKV task CSVs.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path.cwd(),
        help="Directory where canonical ZIKV CSVs will be written.",
    )
    args = parser.parse_args()

    outputs = build_zika_rows()
    for filename, rows in outputs.items():
        destination = args.outdir / filename
        write_rows(destination, rows)
        print(f"Wrote {len(rows)} rows to {destination}")


if __name__ == "__main__":
    main()
