import csv
from decimal import Decimal
from pathlib import Path


def mean(values):
    values = [float(v) for v in values if v is not None and v != ""]
    return sum(values) / len(values)


def exact_mean(values):
    values = [Decimal(str(v)) for v in values if v is not None and v != ""]
    return float(sum(values) / Decimal(len(values)))


def write_rows(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["mutation", "mutated_sequence", "DMS_score"])
        writer.writeheader()
        writer.writerows(rows)


def read_rows(path, strip_terminal_stop=False):
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            mutated_sequence = row["mutated_sequence"]
            if strip_terminal_stop:
                mutated_sequence = mutated_sequence.rstrip("*")
            rows.append(
                {
                    "mutation": row["mutation"],
                    "mutated_sequence": mutated_sequence,
                    "DMS_score": float(row["DMS_score"]),
                }
            )
    return rows


def parse_fasta_sequence(path):
    sequence = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            sequence.append(line)
    return "".join(sequence)


def mutate_sequence(reference_sequence, site, mutant):
    return reference_sequence[: site - 1] + mutant + reference_sequence[site:]


def infer_reference_sequence_from_task_rows(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        first_row = next(reader)

    mutation = first_row["mutation"]
    mutated_sequence = first_row["mutated_sequence"]
    site = int(mutation[1:-1])
    wildtype = mutation[0]
    return mutate_sequence(mutated_sequence, site, wildtype)


TASKS_DIR = Path(__file__).resolve().parent.parent
RAW_CSV_DIR = TASKS_DIR / "raw_csv"


def raw_csv_path(*parts):
    return RAW_CSV_DIR.joinpath(*parts)
