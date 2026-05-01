#!/usr/bin/env python3

import argparse
import csv
from collections import OrderedDict, defaultdict
from pathlib import Path
from statistics import fmean

from processing_cleaned.common import mutate_sequence, raw_csv_path, write_rows


RBD_REGION_START = 331
RBD_REGION_END = 531

RBD_BINDING_OUTPUTS = OrderedDict(
    [
        ("Alpha", "SARS_binding_Alpha_RBD.csv"),
        ("Beta", "SARS_binding_Beta_RBD.csv"),
        ("Delta", "SARS_binding_Delta_RBD.csv"),
        ("Eta", "SARS_binding_Eta_RBD.csv"),
        ("Omicron_BA1", "SARS_binding_Omicron_BA1_RBD.csv"),
        ("Omicron_BA286", "SARS_binding_Omicron_BA286_RBD.csv"),
        ("Omicron_BA2", "SARS_binding_Omicron_BA2_RBD.csv"),
        ("Omicron_BQ11", "SARS_binding_Omicron_BQ11_RBD.csv"),
        ("Omicron_EG5", "SARS_binding_Omicron_EG5_RBD.csv"),
        ("Omicron_FLip", "SARS_binding_Omicron_FLip_RBD.csv"),
        ("Omicron_XBB15", "SARS_binding_Omicron_XBB15_RBD.csv"),
        ("Wuhan-Hu-1", "SARS_binding_Wuhan_Hu_1_RBD.csv"),
    ]
)

RBD_EXPRESSION_OUTPUTS = OrderedDict(
    [
        ("Alpha", "SARS_expression_Alpha_RBD.csv"),
        ("Beta", "SARS_expression_Beta_RBD.csv"),
        ("Delta", "SARS_expression_Delta_RBD.csv"),
        ("Eta", "SARS_expression_Eta_RBD.csv"),
        ("Omicron_BA1", "SARS_expression_Omicron_BA1_RBD.csv"),
        ("Omicron_BA286", "SARS_expression_Omicron_BA286_RBD.csv"),
        ("Omicron_BA2", "SARS_expression_Omicron_BA2_RBD.csv"),
        ("Omicron_BQ11", "SARS_expression_Omicron_BQ11_RBD.csv"),
        ("Omicron_EG5", "SARS_expression_Omicron_EG5_RBD.csv"),
        ("Omicron_FLip", "SARS_expression_Omicron_FLip_RBD.csv"),
        ("Omicron_XBB15", "SARS_expression_Omicron_XBB15_RBD.csv"),
        ("Wuhan-Hu-1", "SARS_expression_Wuhan_Hu_1_RBD.csv"),
    ]
)


def parse_multi_fasta(path):
    sequences = OrderedDict()
    current_name = None
    current_chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(current_chunks)
                current_name = line[1:]
                current_chunks = []
                continue
            current_chunks.append(line)
    if current_name is not None:
        sequences[current_name] = "".join(current_chunks)
    return sequences


def load_csv_rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def patch_reference_sequence(reference_sequence, rows, site_key="site", wildtype_key="wildtype"):
    patched = reference_sequence
    for row in rows:
        site = int(row[site_key])
        wildtype = row[wildtype_key]
        if patched[site - 1] != wildtype:
            patched = patched[: site - 1] + wildtype + patched[site:]
    return patched


def build_rbd_mutated_sequence(reference_sequence, site, mutant):
    return (
        reference_sequence[RBD_REGION_START - 1 : site - 1]
        + mutant
        + reference_sequence[site:RBD_REGION_END]
    )


def build_rbd_rows():
    rows = load_csv_rows(raw_csv_path("SARS-CoV-2", "original", "mutation_effect_RBD.csv"))
    for row in rows:
        row["ACE2 binding"] = row["delta_bind"]
        row["expression"] = row["delta_expr"]
    reference_sequences = parse_multi_fasta(
        raw_csv_path("SARS-CoV-2", "original", "reference.fasta")
    )

    rows_by_target = OrderedDict()
    for row in rows:
        rows_by_target.setdefault(row["target"], []).append(row)

    binding_outputs = {}
    expression_outputs = {}
    for target, target_rows in rows_by_target.items():
        reference_sequence = patch_reference_sequence(reference_sequences[target], target_rows)
        binding_rows = []
        expression_rows = []
        for row in target_rows:
            site = int(row["site"])
            mutation = f"{row['wildtype']}{site - RBD_REGION_START + 1}{row['mutant']}"
            mutated_sequence = build_rbd_mutated_sequence(reference_sequence, site, row["mutant"])

            if row["ACE2 binding"] != "":
                binding_rows.append(
                    {
                        "mutation": mutation,
                        "mutated_sequence": mutated_sequence,
                        "DMS_score": float(row["ACE2 binding"]),
                    }
                )
            if row["expression"] != "":
                expression_rows.append(
                    {
                        "mutation": mutation,
                        "mutated_sequence": mutated_sequence,
                        "DMS_score": float(row["expression"]),
                    }
                )

        binding_outputs[RBD_BINDING_OUTPUTS[target]] = binding_rows
        expression_outputs[RBD_EXPRESSION_OUTPUTS[target]] = expression_rows
    return binding_outputs | expression_outputs


def build_xbb15_rbd_rows(reference_sequences):
    rows = load_csv_rows(raw_csv_path("SARS-CoV-2", "original", "XBB15_RBD.csv"))
    reference_sequence = patch_reference_sequence(reference_sequences["Omicron_XBB15"], rows)

    outputs = {
        "SARS_sera_escape_Omicron_XBB15_RBD.csv": [],
        "SARS_cell_entry_Omicron_XBB15_RBD.csv": [],
    }
    for row in rows:
        site = int(row["site"])
        mutation = f"{row['wildtype']}{site - RBD_REGION_START + 1}{row['mutant']}"
        mutated_sequence = build_rbd_mutated_sequence(reference_sequence, site, row["mutant"])

        if row["human sera escape"] != "":
            outputs["SARS_sera_escape_Omicron_XBB15_RBD.csv"].append(
                {
                    "mutation": mutation,
                    "mutated_sequence": mutated_sequence,
                    "DMS_score": float(row["human sera escape"]),
                }
            )
        if row["spike mediated entry"] != "":
            outputs["SARS_cell_entry_Omicron_XBB15_RBD.csv"].append(
                {
                    "mutation": mutation,
                    "mutated_sequence": mutated_sequence,
                    "DMS_score": float(row["spike mediated entry"]),
                }
            )
    return outputs


def build_full_spike_rows(source_name, target_name, binding_output_name, cell_entry_output_name):
    rows = load_csv_rows(raw_csv_path("SARS-CoV-2", "original", source_name))
    reference_sequences = parse_multi_fasta(
        raw_csv_path("SARS-CoV-2", "original", "reference.fasta")
    )
    reference_sequence = patch_reference_sequence(reference_sequences[target_name], rows)

    outputs = {
        binding_output_name: [],
        cell_entry_output_name: [],
    }
    if source_name == "XBB15_full.csv":
        outputs["SARS_sera_escape_Omicron_XBB15.csv"] = []

    for row in rows:
        site = int(row["site"])
        mutation = f"{row['wildtype']}{site}{row['mutant']}"
        mutated_sequence = mutate_sequence(reference_sequence, site, row["mutant"])

        if row["ACE2 binding"] != "":
            outputs[binding_output_name].append(
                {
                    "mutation": mutation,
                    "mutated_sequence": mutated_sequence,
                    "DMS_score": float(row["ACE2 binding"]),
                }
            )
        if row["spike mediated entry"] != "":
            outputs[cell_entry_output_name].append(
                {
                    "mutation": mutation,
                    "mutated_sequence": mutated_sequence,
                    "DMS_score": float(row["spike mediated entry"]),
                }
            )
        if source_name == "XBB15_full.csv" and row["human sera escape"] != "":
            outputs["SARS_sera_escape_Omicron_XBB15.csv"].append(
                {
                    "mutation": mutation,
                    "mutated_sequence": mutated_sequence,
                    "DMS_score": float(row["human sera escape"]),
                }
            )
    return outputs


def build_wuhan_antibody_escape_rows(reference_sequences):
    rows = load_csv_rows(raw_csv_path("SARS-CoV-2", "mAb", "high_antibody_results.csv"))
    reference_sequence = reference_sequences["Wuhan-Hu-1"]
    grouped_scores = defaultdict(list)

    for row in rows:
        site = int(row["site"])
        mutation = f"{row['wildtype']}{site - RBD_REGION_START + 1}{row['mutation']}"
        mutated_sequence = build_rbd_mutated_sequence(reference_sequence, site, row["mutation"])
        key = (mutation, mutated_sequence)
        grouped_scores[key].append(float(row["escape_score"]))

    output_rows = []
    for (mutation, mutated_sequence), scores in grouped_scores.items():
        output_rows.append(
            {
                "mutation": mutation,
                "mutated_sequence": mutated_sequence,
                "DMS_score": fmean(scores),
            }
        )
    output_rows.sort(key=lambda row: (row["mutation"], row["mutated_sequence"]))
    return {"SARS_antibody_escape_Wuhan_Hu_1.csv": output_rows}


def build_sars_rows():
    reference_sequences = parse_multi_fasta(
        raw_csv_path("SARS-CoV-2", "original", "reference.fasta")
    )

    outputs = OrderedDict()
    outputs.update(build_rbd_rows())
    outputs.update(build_xbb15_rbd_rows(reference_sequences))
    outputs.update(
        build_full_spike_rows(
            "XBB15_full.csv",
            "Omicron_XBB15",
            "SARS_binding_Omicron_XBB15.csv",
            "SARS_cell_entry_Omicron_XBB15.csv",
        )
    )
    outputs.update(
        build_full_spike_rows(
            "BA2.csv",
            "Omicron_BA2",
            "SARS_binding_Omicron_BA2.csv",
            "SARS_cell_entry_Omicron_BA2.csv",
        )
    )
    outputs.update(build_wuhan_antibody_escape_rows(reference_sequences))
    return outputs


def main():
    parser = argparse.ArgumentParser(
        description="Export canonical SARS-CoV-2 task CSVs from checked-in source tables."
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path.cwd(),
        help="Directory where canonical SARS-CoV-2 CSVs will be written.",
    )
    args = parser.parse_args()

    outputs = build_sars_rows()
    for filename, rows in outputs.items():
        destination = args.outdir / filename
        write_rows(destination, rows)
        print(f"Wrote {len(rows)} rows to {destination}")


if __name__ == "__main__":
    main()
