Cleaned task exporters live here.

This package is intentionally flat: exporter modules live directly under `processing_cleaned/` and are discovered through the package registry.

## Run

```bash
python -m processing_cleaned list
python -m processing_cleaned export --outdir ./exports
```

`export_all.py` provides the same interface from the repository root.

## Structure

- `common.py`: shared parsing, mutation, averaging, and CSV-writing helpers
- `registry.py`: exporter registry used by the package CLI
- `export_*.py`: runnable exporter modules

## Implemented exporters

- `export_hbv_fitness.py`
- `export_hiv_sera_escape_b520.py`
- `export_hiv_mab_tasks.py`
- `export_rabv_tasks.py`
- `export_nipah_tasks.py`
- `export_lassa_tasks.py`
- `export_zika_tasks.py`
- `export_h3n2_mc22_tasks.py`
- `export_h3n2_hk19_tasks.py`
- `export_h5n1_tasks.py`
- `export_proteingym_tasks.py`
- `export_sars_task_tables.py`
- `export_sars_kp311_tasks.py`

## Conventions

- Exporters read frozen inputs from `raw_csv/`.
- Outputs use the canonical schema `mutation`, `mutated_sequence`, `DMS_score`.
- Legacy benchmark behavior is preserved where reproducibility matters more than semantic cleanup.