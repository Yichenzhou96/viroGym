#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

from processing_cleaned.common import raw_csv_path, write_rows


PROTEINGYM_EXPORTS = {
    "RDRP_I33A0_Li_2023.csv": "RDRP_I33A0_Li_2023.csv",
    "A0A2Z5U3Z0_9INFA_Wu_2014.csv": "A0A2Z5U3Z0_9INFA_Wu_2014.csv",
    "A0A2Z5U3Z0_9INFA_Doud_2016.csv": "FLU_viral_growth_H1N1.csv",
    "A4D664_9INFA_Soh_2019.csv": "FLU_viral_growth_H2N1.csv",
    "CAPSD_AAV2S_Sinai_2021.csv": "AAV2_viral_growth.csv",
    "C6KNH7_9INFA_Lee_2018.csv": "FLU_viral_growth_H3N2.csv",
    "Q2N0S5_9HIV1_Haddox_2018.csv": "HIV_viral_growth.csv",
    "A0A192B1T2_9HIV1_Haddox_2018.csv": "HIV_viral_growth_BF520.csv",
    "ENV_HV1BR_Haddox_2016.csv": "HIV_viral_growth_BRU_LAI.csv",
    "REV_HV1H2_Fernandes_2016.csv": "HIV_viral_growth_HXB2.csv",
    "ENV_HV1B9_DuenasDecamp_2016.csv": "HIV_viral_growth_strain896.csv",
    "I6TAH8_I68A0_Doud_2015.csv": "I6TAH8_I68A0_Doud_2015.csv",
    "POLG_CXB3N_Mattenberger_2021.csv": "CXB3_viral_growth.csv",
    "POLG_DEN26_Suphatrakul_2023.csv": "DEN_viral_growth.csv",
    "POLG_HCVJF_Qi_2014.csv": "HCV_viral_growth_JFH_1.csv",
    "POLG_PESV_Tsuboyama_2023_2MXD.csv": "PESV_stability.csv",
    "NCAP_I34A1_Doud_2015.csv": "NCAP_I34A1_Doud_2015.csv",
    "NRAM_I33A0_Jiang_2016.csv": "NRAM_I33A0_Jiang_2016.csv",
    "PA_I34A1_Wu_2015.csv": "PA_I34A1_Wu_2015.csv",
    "R1AB_SARS2_Flynn_2022.csv": "SARS_viral_growth_Wuhan_Hu_1.csv",
    "TAT_HV1BR_Fernandes_2016.csv": "TAT_HV1BR_Fernandes_2016.csv",
}


def convert_rows(reader):
    rows = []
    for row in reader:
        rows.append(
            {
                "mutation": row["mutant"],
                "mutated_sequence": row["mutated_sequence"],
                "DMS_score": float(row["DMS_score"]),
            }
        )
    return rows


def load_local_rows(path):
    with open(path, newline="") as handle:
        return convert_rows(csv.DictReader(handle))

def build_proteingym_rows():
    outputs = {}
    for source_name, output_name in PROTEINGYM_EXPORTS.items():
        outputs[output_name] = load_local_rows(raw_csv_path("ProteinGym", source_name))
    return outputs


def main():
    parser = argparse.ArgumentParser(
        description="Export canonical benchmark CSVs sourced from local ProteinGym data."
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path.cwd(),
        help="Directory where canonical ProteinGym-backed CSVs will be written.",
    )
    args = parser.parse_args()

    outputs = build_proteingym_rows()
    for filename, rows in outputs.items():
        destination = args.outdir / filename
        write_rows(destination, rows)
        print(f"Wrote {len(rows)} rows to {destination}")


if __name__ == "__main__":
    main()
