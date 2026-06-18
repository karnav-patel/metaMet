import os
import csv
import codecs
import re
import html
import sys
from pathlib import Path
from collections import defaultdict

def get_project_root():
    current_file_path = Path(__file__).resolve()
    for parent in current_file_path.parents:
        if parent.name == "metaMet":
            return parent
    raise RuntimeError("Project root folder 'metaMet' not found.")

# locate and import your config
project_root = get_project_root()
data_preproc_dir = project_root / "data_preprocessing"
sys.path.insert(0, str(data_preproc_dir))

import config.config as config
import subprocess

def check_path_exists(path, description):
    """
    Prints and returns whether a given file path exists.
    """
    p = Path(path)
    if p.exists():
        print(f"[OK] {description} exists: {p}")
        return True
    else:
        print(f"[ERROR] {description} does not exist: {p}", file=sys.stderr)
        return False

def check_all_paths(paths_with_desc):
    """
    Given a list of (description, path) tuples, checks each and
    returns True only if they all exist.
    """
    all_good = True
    for desc, p in paths_with_desc:
        if not check_path_exists(p, desc):
            all_good = False
    return all_good

def run_script(script_path, description):
    """
    Runs a Python script and prints status. Returns True on success.
    """
    print(f"[RUN] {description}: {script_path}")
    try:
        subprocess.run([sys.executable, str(script_path)], check=True)
        print(f"[OK] Completed: {description}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] {description} failed with exit code {e.returncode}", file=sys.stderr)
        return False

def main():
    # 1) Pre-check inputs for the prepare step
    prepare_checks = [
        ("BRENDA kcat CSV", config.brenda_kcat_csv),
        ("MetaCyc kcat CSV", config.metacyc_kcat_csv),
        ("SABIO kcat TSV", config.sabio_kcat_csv),
        ("Master EC list", config.merged_extracted_ids_output_path),
        ("MetaCyc reactions-with-SMILES", config.metacyc_reaction_csv_with_smiles),
        ("BRENDA reactions CSV", config.brenda_reaction_csv),
        ("KEGG reactions CSV", config.kegg_reaction_csv),
        ("MetaCyc reactions CSV", config.metacyc_reaction_csv),
    ]
    if not check_all_paths(prepare_checks):
        sys.exit(1)

    # RUN FIRST: fetch EC → protein sequences
    # if not check_path_exists(config.predictkcat_protein_sequences, "Fetch EC protein sequences script"):
    #     sys.exit(1)
    # if not run_script(config.predictkcat_protein_sequences, "Fetch EC protein sequences"):
    #     sys.exit(1)

    # Ensure the sequence CSV produced by the fetch step exists
    if not check_path_exists(config.predictkcat_with_ec_protein_sequences, "EC→protein sequence CSV"):
        sys.exit(1)


    # 2) Run the prepare script
    if not run_script(config.predictkcat_prepare, "Prepare kcat inputs"):
        sys.exit(1)

    # 3) Pre-check inputs for the predict step
    predict_checks = [
        ("Predict input TSV", config.predictkcat_input_tsv),
        ("Predict input-with-EC CSV", config.predictkcat_input_with_ec),
        ("Predict failed ECs TXT", config.predictkcat_failed_kcat_ids_output_path),
    ]
    if not check_all_paths(predict_checks):
        sys.exit(1)

    # 4) Run the predict script
    if not run_script(config.predictkcat_predict, "Run kcat prediction"):
        sys.exit(1)

    # 5) Final output verification
    if check_path_exists(config.predictkcat_output_tsv, "Predict output TSV"):
        print("[OK] All steps completed successfully.")
        sys.exit(0)
    else:
        print("[ERROR] Final output missing; prediction failed.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
