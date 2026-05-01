"""Exporter registry for cleaned task-processing scripts."""

from dataclasses import dataclass
from importlib import import_module
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class ExporterSpec:
    name: str
    module: str
    build_attr: str
    description: str
    output_name: Optional[str] = None


EXPORTERS = (
    ExporterSpec(
        name="hbv_fitness",
        module="processing_cleaned.export_hbv_fitness",
        build_attr="build_hbv_rows",
        output_name="HBV_fitness.csv",
        description="HBV fitness export.",
    ),
    ExporterSpec(
        name="hiv_sera_escape_b520",
        module="processing_cleaned.export_hiv_sera_escape_b520",
        build_attr="build_hiv_b520_rows",
        output_name="HIV_sera_escape_B520.csv",
        description="HIV BF520 serum-escape export.",
    ),
    ExporterSpec(
        name="hiv_mab",
        module="processing_cleaned.export_hiv_mab_tasks",
        build_attr="build_hiv_mab_rows",
        description="HIV mAb antibody-escape and cell-entry exports.",
    ),
    ExporterSpec(
        name="rabv",
        module="processing_cleaned.export_rabv_tasks",
        build_attr="build_rabv_rows",
        description="RABV antibody-escape and cell-entry exports.",
    ),
    ExporterSpec(
        name="nipah",
        module="processing_cleaned.export_nipah_tasks",
        build_attr="build_nipah_rows",
        description="Nipah antibody-escape, binding, and cell-entry exports.",
    ),
    ExporterSpec(
        name="lassa",
        module="processing_cleaned.export_lassa_tasks",
        build_attr="build_lassa_rows",
        description="LASV antibody-escape and cell-entry exports.",
    ),
    ExporterSpec(
        name="influenza_h5n1",
        module="processing_cleaned.export_h5n1_tasks",
        build_attr="build_h5n1_rows",
        description="Influenza H5N1 cell-entry, sera-escape, and stability exports.",
    ),
    ExporterSpec(
        name="influenza_h3n2_mc22",
        module="processing_cleaned.export_h3n2_mc22_tasks",
        build_attr="build_h3n2_mc22_rows",
        description="Influenza H3N2 MC22 exports.",
    ),
    ExporterSpec(
        name="influenza_h3n2_hk19",
        module="processing_cleaned.export_h3n2_hk19_tasks",
        build_attr="build_hk19_rows",
        description="Influenza H3N2 HK19 exports.",
    ),
    ExporterSpec(
        name="proteingym",
        module="processing_cleaned.export_proteingym_tasks",
        build_attr="build_proteingym_rows",
        description="ProteinGym-backed benchmark exports.",
    ),
    ExporterSpec(
        name="sars_cov_2",
        module="processing_cleaned.export_sars_task_tables",
        build_attr="build_sars_rows",
        description="Canonical SARS-CoV-2 exports.",
    ),
    ExporterSpec(
        name="sars_cov_2_kp311",
        module="processing_cleaned.export_sars_kp311_tasks",
        build_attr="build_kp311_rows",
        description="Canonical KP.3.1.1 SARS-CoV-2 exports.",
    ),
    ExporterSpec(
        name="zika",
        module="processing_cleaned.export_zika_tasks",
        build_attr="build_zika_rows",
        description="ZIKV antibody-escape and viral-growth exports.",
    ),
)


def iter_exporters():
    return EXPORTERS


def build_export_rows(spec):
    module = import_module(spec.module)
    rows = getattr(module, spec.build_attr)()
    if spec.output_name is not None:
        return {spec.output_name: rows}
    return rows


def export_to_directory(spec, outdir, write_rows):
    outdir = Path(outdir)
    outputs = build_export_rows(spec)
    written = []
    for filename, rows in outputs.items():
        destination = outdir / filename
        write_rows(destination, rows)
        written.append(destination)
    return written
