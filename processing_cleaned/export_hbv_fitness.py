#!/usr/bin/env python3

import argparse
import csv
from collections import OrderedDict
from pathlib import Path

from processing_cleaned.common import mean, raw_csv_path, write_rows


def build_hbv_rows():
    input_path = raw_csv_path("HBV", "DMS_HBV_DNA.csv")

    rows = []
    max_site = 0
    with open(input_path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            site = int(float(row["site"]))
            max_site = max(max_site, site)
            rows.append(
                {
                    "site": site,
                    "aa": row["AA"],
                    "wt_aa": row["wtAA"],
                    "score": row["logmfactor"],
                }
            )

    ref_index = [None] * max_site
    for idx, row in enumerate(rows):
        ref_index[row["site"] - 1] = idx

    ref_sequence = "".join(rows[idx]["wt_aa"] for idx in ref_index)

    aggregated = OrderedDict()
    for row in rows:
        mutation = f"{row['wt_aa']}{row['site']}{row['aa']}"
        mutated_sequence = ref_sequence[: row["site"] - 1] + row["aa"] + ref_sequence[row["site"] :]
        key = (mutation, mutated_sequence)
        aggregated.setdefault(key, []).append(float(row["score"]))

    return [
        {
            "mutation": mutation,
            "mutated_sequence": mutated_sequence,
            "DMS_score": mean(scores),
        }
        for (mutation, mutated_sequence), scores in aggregated.items()
    ]


def main():
    parser = argparse.ArgumentParser(description="Export canonical HBV_fitness.csv rows.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path.cwd() / "HBV_fitness.csv",
        help="Destination CSV path.",
    )
    args = parser.parse_args()

    rows = build_hbv_rows()
    write_rows(args.output, rows)
    print(f"Wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
