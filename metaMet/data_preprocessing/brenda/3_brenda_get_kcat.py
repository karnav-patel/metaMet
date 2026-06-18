#!/usr/bin/env python3
import sys
import re
from pathlib import Path

import pandas as pd
import ijson

# -----------------------------------------------------------------------------
# Helper functions to parse pH and temperature from comment text
# -----------------------------------------------------------------------------
def parse_ph(comment):
    """
    Look for 'pH <number>' in the comment (e.g. "pH 10.0 <49>").
    Returns a float or None if not found.
    """
    if not comment:
        return None
    match = re.search(r"pH\s*([\d]+(?:\.[\d]+)?)", comment)
    return float(match.group(1)) if match else None


def parse_temperature(comment):
    """
    Look for '<number> °C' in the comment (e.g. "at 65°C").
    Returns a float or None if not found.
    """
    if not comment:
        return None
    match = re.search(r"(\d+(?:\.[\d]+)?)\s*°C", comment)
    return float(match.group(1)) if match else None


# -----------------------------------------------------------------------------
# Dynamically locate project root and load configuration
# -----------------------------------------------------------------------------
def get_project_root():
    current_file = Path(__file__).resolve()
    for parent in current_file.parents:
        if parent.name == "metaMet":
            return parent
    raise RuntimeError("Project root folder 'metaMet' not found.")

project_root = get_project_root()
data_preproc_dir = project_root / "data_preprocessing"
sys.path.insert(0, str(data_preproc_dir))

# Import paths from config/config.py
import config.config as config


def extract_brenda_kcat_data():
    # Paths from config
    brenda_json = config.brenda_json_path
    extracted_ec_file = config.brenda_extracted_ids_output_path
    output_csv = config.brenda_kcat_csv

    # Read EC numbers
    with open(extracted_ec_file, 'r') as f:
        ec_numbers = [line.strip() for line in f if line.strip()]

    records = []
    for ec in ec_numbers:
        print(f"Processing EC number: {ec}")
        with open(brenda_json, 'r') as json_file:
            ec_path = f"data.{ec}"
            try:
                for item in ijson.items(json_file, ec_path):
                    if 'turnover_number' not in item:
                        continue
                    for tn in item['turnover_number']:
                        value = tn.get('num_value')
                        if value is None:
                            continue

                        # Core fields
                        substrate = tn.get('value')
                        comment = tn.get('comment', '')
                        refs = tn.get('references', [])
                        orgs = tn.get('organisms', [])
                        prots = tn.get('proteins', [])

                        kcat = float(value)
                        ph_val = parse_ph(comment)
                        temp_val = parse_temperature(comment)

                        records.append({
                            'ec_number': ec,
                            'substrate': substrate,
                            'kcat': kcat,
                            'comment': comment,
                            'pH': ph_val,
                            'temperature_C': temp_val,
                            'references': ';'.join(refs),
                            'organisms': ';'.join(orgs),
                            'proteins': ';'.join(prots),
                        })
            except ijson.JSONError as e:
                print(f"JSON error for {ec}: {e}")
            except StopIteration:
                print(f"No data for EC number {ec}")

    # Write to CSV
    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    print(f"Kcat data written to {output_csv}")


if __name__ == '__main__':
    extract_brenda_kcat_data()
