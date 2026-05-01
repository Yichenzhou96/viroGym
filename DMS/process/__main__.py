#!/usr/bin/env python3

import argparse
from pathlib import Path

from processing_cleaned.common import write_rows
from processing_cleaned.registry import export_to_directory, iter_exporters


def build_parser():
    parser = argparse.ArgumentParser(description="Run cleaned task exporters.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    list_parser = subparsers.add_parser("list", help="List available exporters.")
    list_parser.set_defaults(func=run_list)

    export_parser = subparsers.add_parser("export", help="Run one or more exporters.")
    export_parser.add_argument(
        "--outdir",
        type=Path,
        default=Path.cwd(),
        help="Directory where exported CSV files will be written.",
    )
    export_parser.add_argument(
        "names",
        nargs="*",
        help="Exporter names to run. Omit to run all exporters.",
    )
    export_parser.set_defaults(func=run_export)
    return parser


def run_list(_args):
    for spec in iter_exporters():
        print(f"{spec.name}\t{spec.description}")


def run_export(args):
    selected = set(args.names)
    specs = [spec for spec in iter_exporters() if not selected or spec.name in selected]
    if selected:
        found = {spec.name for spec in specs}
        missing = sorted(selected - found)
        if missing:
            raise SystemExit(f"Unknown exporters: {', '.join(missing)}")

    for spec in specs:
        for destination in export_to_directory(spec, args.outdir, write_rows):
            print(f"{spec.name}\t{destination}")


def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
