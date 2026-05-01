#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

from processing_cleaned.common import mutate_sequence, raw_csv_path, write_rows


H5N1_REFERENCE_SEQUENCE = (
    "MENIVLLLAIVSLVKSDQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKTHNGKLCDLNGVKPLILKDCSVAGWLLGNPMCDEFIRVPEWSYIVERANPANDLCYPGSLNDYEELKHMLSRINHFEKILIIPKSSWPNHETSLGVSAACPYQGAPSFFRNVVWLIKKNDAYPTIKISYNNTNREDLLILWGIHHSNNAEEQTNLYKNPTTYISVGTSTLNQRLAPKIATRSQVNGQRGRMDFFWTILKPDDAIHFESNGNFIAPEYAYKIVKKGDSTIMKSGVEYGHCNTKCQTPVGAINSSMPFHNIHPLTIGECPKYVKSNKLVLATGLRNSPLREKRRKRGLFGAIAGFIEGGWQGMVDGWYGYHHSNEQGSGYAADKESTQKAIDGVTNKVNSIIDKMNTQFEAVGREFNNLERRIENLNKKMEDGFLDVWTYNAELLVLMENERTLDFHDSNVKNLYDKVRLQLRDNAKELGNGCFEFYHKCDNECMESVRNGTYDYPQYSEEARLKREEISGVKLESVGTYQILSIYSTAASSLALAIMMAGLSLWMCSNGSLQCRICI"
)


def build_h5n1_rows():
    source_path = raw_csv_path("influenza", "phenotypes.csv")

    outputs = {
        "FLU_sera_escape_H5N1.csv": [],
        "FLU_cell_entry_H5N1.csv": [],
        "FLU_stability_H5N1.csv": [],
    }
    score_columns = {
        "FLU_sera_escape_H5N1.csv": "species sera escape",
        "FLU_cell_entry_H5N1.csv": "entry in 293T cells",
        "FLU_stability_H5N1.csv": "stability",
    }

    with open(source_path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            site = int(row["sequential_site"])
            mutation = f"{row['wildtype']}{site}{row['mutant']}"
            mutated_sequence = mutate_sequence(H5N1_REFERENCE_SEQUENCE, site, row["mutant"])

            for output_name, column_name in score_columns.items():
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
    parser = argparse.ArgumentParser(description="Export canonical H5N1 influenza task CSVs.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path.cwd(),
        help="Directory where canonical H5N1 influenza CSVs will be written.",
    )
    args = parser.parse_args()

    outputs = build_h5n1_rows()
    for filename, rows in outputs.items():
        destination = args.outdir / filename
        write_rows(destination, rows)
        print(f"Wrote {len(rows)} rows to {destination}")


if __name__ == "__main__":
    main()
