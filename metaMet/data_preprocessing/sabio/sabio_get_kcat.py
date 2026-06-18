import os
import re

import sys
from pathlib import Path

def get_project_root():
    current_file_path = Path(__file__).resolve()
    for parent in current_file_path.parents:
        if parent.name == "metaMet":
            return parent
    raise RuntimeError("Project root folder 'metaMet' not found.")

project_root = get_project_root()
data_preproc_dir = project_root / "data_preprocessing"
sys.path.insert(0, str(data_preproc_dir))

import config.config as config 

"""
Combined pipeline for:
1. Extracting EC numbers
2. Fetching SABIO kinetics TSV data per EC
3. Generating unisubstrate Kcat/Km TSV
4. Cleaning and selecting max Kcat values

Author: Adapted to use project config
Date: 2025-07-07
"""
import csv
import requests

# --- Paths from config.py ---
EC_LIST_FILE = config.merged_extracted_ids_output_path  # Path to merged EC numbers file
SABIO_EC_DIR = config.sabio_temp_data_ec_number  # Subfolder for per-EC text files
UNISUBSTRATE_TSV = config.sabio_temp_data_kcat_unisubstrate  # TSV for Kcat/Km pairs
CLEAN_TSV = config.sabio_kcat_csv                       # Final cleaned Kcat TSV

if not EC_LIST_FILE.exists():
    raise FileNotFoundError(f"EC list not found: {EC_LIST_FILE}")

# SABIO REST endpoint
SABIO_QUERY_URL = (
    'http://sabiork.h-its.org/sabioRestWebServices/kineticlawsExportTsv'
)


def get_ec_list():
    """
    Read EC numbers from the merged EC list file, one per line.
    """
    if not EC_LIST_FILE.exists():
        raise FileNotFoundError(f"EC list not found: {EC_LIST_FILE}")
    with EC_LIST_FILE.open('r') as infile:
        ec_list = [line.strip() for line in infile if line.strip()]
    print(f"Loaded {len(ec_list)} EC numbers from {EC_LIST_FILE}")
    return ec_list


def fetch_sabio_info(ec_list):
    """
    For each EC number, request SABIO kinetics and write to individual files.
    """
    SABIO_EC_DIR.mkdir(parents=True, exist_ok=True)
    for idx, ec_num in enumerate(ec_list, start=1):
        print(f"Fetching {idx}/{len(ec_list)}: EC {ec_num}")
        params = {
            'fields[]': [
                'EntryID', 'Substrate', 'EnzymeType', 'PubMedID',
                'Organism', 'UniprotID', 'ECNumber', 'Parameter'
            ],
            'q': f"ECNumber:{ec_num}"
        }
        resp = requests.post(SABIO_QUERY_URL, params=params)
        if resp.ok and resp.text:
            out_file = SABIO_EC_DIR / f"{ec_num}.txt"
            out_file.write_text(resp.text)
        else:
            print(f"Warning: no data for EC {ec_num}")


def generate_unisubstrate_tsv():
    """
    Parse each EC file and write Kcat/Km pairs for unisubstrate kinetics.
    """
    # Ensure parent directory exists
    UNISUBSTRATE_TSV.parent.mkdir(parents=True, exist_ok=True)
    with UNISUBSTRATE_TSV.open('w', newline='') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        writer.writerow([
            'EntryID', 'Type', 'ECNumber', 'Substrate', 'EnzymeType',
            'PubMedID', 'Organism', 'UniprotID', 'Value', 'Unit'
        ])
        for txt in SABIO_EC_DIR.iterdir():
            if txt.suffix.lower() != '.txt':
                continue
            lines = txt.read_text(encoding='utf-8').splitlines()
            # Build Km lookup per EntryID
            km_lookup = {}
            for row in lines[1:]:
                cols = row.split('\t')
                if len(cols) > 7 and cols[7] == 'Km':
                    key = cols[0]
                    km_lookup.setdefault(key, []).append((cols[9], cols[-1]))
            # Pair kcat with Km
            for row in lines[1:]:
                cols = row.split('\t')
                try:
                    if cols[7] == 'kcat' and len(cols) > 10 and cols[10]:
                        entry = cols[0]
                        for km_val, km_unit in km_lookup.get(entry, []):
                            writer.writerow([
                                entry, cols[7], cols[6], km_val,
                                cols[2], cols[3], cols[4], cols[5],
                                cols[10], cols[-1]
                            ])
                except Exception:
                    continue
    print(f"Unisubstrate TSV generated at {UNISUBSTRATE_TSV}")


def clean_kcat_data():
    """
    From the unisubstrate TSV, select the maximum Kcat per unique entry and
    unify units to 's^(-1)'.
    """
    if not UNISUBSTRATE_TSV.exists():
        raise FileNotFoundError(f"Input TSV not found: {UNISUBSTRATE_TSV}")
    # Read all records
    with UNISUBSTRATE_TSV.open('r', encoding='utf-8') as in_f:
        reader = csv.reader(in_f, delimiter='\t')
        header = next(reader)
        records = list(reader)

    # Determine max Kcat per unique key
    unique = {}
    for cols in records:
        key = tuple(cols[:7])  # first seven columns
        try:
            val = float(cols[8])
        except ValueError:
            continue
        unit = cols[9]
        if key not in unique or val > unique[key][0]:
            unique[key] = (val, unit)

    # Write cleaned output
    CLEAN_TSV.parent.mkdir(parents=True, exist_ok=True)
    with CLEAN_TSV.open('w', newline='') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        writer.writerow(header)
        for key, (val, unit) in unique.items():
            # normalize unit
            if unit in ['mol*s^(-1)*mol^(-1)', 's^(-', '-']:
                unit = 's^(-1)'
            if unit == 's^(-1)':
                writer.writerow(list(key) + [str(val), unit])
    print(f"Cleaned Kcat data written to {CLEAN_TSV}")

def main():
    ecs = get_ec_list()
    fetch_sabio_info(ecs)
    generate_unisubstrate_tsv()
    clean_kcat_data()

if __name__ == '__main__':
    main()
