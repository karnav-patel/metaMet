#!/usr/bin/env python3
import sys
from pathlib import Path

# Dynamically locate project root

def get_project_root():
    current = Path(__file__).resolve()
    for p in current.parents:
        if p.name == "metaMet":
            return p
    raise RuntimeError("Project root folder 'metaMet' not found.")

project_root = get_project_root()
# Ensure data_preprocessing and config are importable
sys.path.insert(0, str(project_root / "data_preprocessing"))
import config.config as config


def load_ecs(path):
    """Read EC numbers from a text file, validating numeric format."""
    ecs = set()
    try:
        with open(path, 'r') as f:
            for line in f:
                ec = line.strip()
                if not ec:
                    continue
                parts = ec.split('.')
                # Only accept if all parts are integers
                if all(p.isdigit() for p in parts):
                    ecs.add(ec)
                else:
                    print(f"Skipping invalid EC format: {ec}")
    except FileNotFoundError:
        print(f"WARNING: Input file not found: {path}")
    return ecs


def ec_sort_key(ec):
    """Convert an EC string like '1.2.3.4' into a tuple of ints for sorting."""
    return tuple(int(n) for n in ec.split('.'))


def main():
    # Gather EC lists from all processing steps
    inputs = [
        config.brenda_extracted_ids_output_path,
        config.kegg_extracted_ids_output_path,
        config.metacyc_extracted_ids_output_path,
    ]
    union_set = set()

    for p in inputs:
        print(f"Loading EC numbers from: {p}")
        ecs = load_ecs(p)
        print(f"  Found {len(ecs)} valid EC numbers.")
        union_set.update(ecs)

    # Sort and write the union
    sorted_ecs = sorted(union_set, key=ec_sort_key)
    out_path = Path(config.merged_extracted_ids_output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as out_f:
        for ec in sorted_ecs:
            out_f.write(ec + "\n")

    print(f"Wrote {len(sorted_ecs)} unique EC numbers to: {out_path}")

if __name__ == '__main__':
    main()
