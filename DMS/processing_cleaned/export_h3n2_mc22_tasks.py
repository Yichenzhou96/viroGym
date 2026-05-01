#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

from processing_cleaned.common import mutate_sequence, raw_csv_path, write_rows


H3N2_MC22_REFERENCE_SEQUENCE = (
    "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGKICNSPHQILDGGNCTLIDALLGDPQCDGFQNKEWDLFVERSRANSSCYPYDVPDYASLRSLVASSGTLEFKNESFNWTGVKQNGTSSACKRGSSSSFFSRLNWLTSLNNIYPAQNVTMPNKEQFDKLYIWGVHHPDTDKNQFSLFAQSSGRITVSTKRSQQAVIPNIGSRPRVRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRITYGACPRYVKQSTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGRGQAADLKSTQAAIDQISGKLNRLIGKTNEKFHQIEKEFSEVEGRVQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNETYDHNVYRDEALNNRFQIKG"
)


OUTPUT_TO_COLUMN = {
    "FLU_cell_entry_H3N2_MC22.csv": "MDCKSIAT1 cell entry",
    "FLU_sera_escape_H3N2_MC22.csv": "sera escape",
    "FLU_stability_H3N2_MC22.csv": "pH stability",
}


def build_h3n2_mc22_rows():
    source_path = raw_csv_path("influenza", "H3N2_MC22", "H3N2.csv")
    outputs = {name: [] for name in OUTPUT_TO_COLUMN}

    with open(source_path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            site = int(row["sequential_site"])
            mutation = f"{row['wildtype']}{site}{row['mutant']}"
            mutated_sequence = mutate_sequence(H3N2_MC22_REFERENCE_SEQUENCE, site, row["mutant"])

            for output_name, column_name in OUTPUT_TO_COLUMN.items():
                score = row[column_name]
                if score == "":
                    continue
                outputs[output_name].append(
                    {
                        "mutation": mutation,
                        "mutated_sequence": mutated_sequence,
                        "DMS_score": float(score),
                    }
                )
    return outputs


def main():
    parser = argparse.ArgumentParser(description="Export canonical H3N2 MC22 influenza task CSVs.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path.cwd(),
        help="Directory where canonical H3N2 MC22 CSVs will be written.",
    )
    args = parser.parse_args()

    outputs = build_h3n2_mc22_rows()
    for filename, rows in outputs.items():
        destination = args.outdir / filename
        write_rows(destination, rows)
        print(f"Wrote {len(rows)} rows to {destination}")


if __name__ == "__main__":
    main()
