#!/usr/bin/env python3
import argparse
import importlib.util
import sys
from pathlib import Path

def load_config_module(config_path):
    """
    Loads a Python config file from the given path as a module.
    """
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    spec = importlib.util.spec_from_file_location("config", str(config_path))
    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)
    return config

def check_path_exists(path, description):
    """
    Checks if a given path (file or directory) exists.
    Returns True if exists, False otherwise.
    """
    p = Path(path)
    if p.exists():
        print(f"[OK] {description} exists: {p}")
        return True
    else:
        print(f"[ERROR] {description} does not exist: {p}")
        return False

def check_parent_exists(path, description):
    """
    Checks if the parent directory of a given file path exists.
    """
    p = Path(path).parent
    if p.exists():
        print(f"[OK] Parent directory for {description} exists: {p}")
        return True
    else:
        print(f"[ERROR] Parent directory for {description} does not exist: {p}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Precheck script to verify configuration paths."
    )
    parser.add_argument(
        "--config", "-c",
        required=True,
        help="Path to the configuration file (Python module) to check."
    )
    args = parser.parse_args()

    # Load the configuration module.
    try:
        config = load_config_module(args.config)
    except Exception as e:
        print("Error loading config module:")
        print(e)
        sys.exit(1)

    # Keep track of whether all checks passed.
    all_ok = True

    # Check that the project root exists.
    if hasattr(config, "project_root"):
        all_ok &= check_path_exists(config.project_root, "Project root")
    else:
        print("[ERROR] Config does not define 'project_root'.")
        all_ok = False

    # Check that the raw data directories exist.
    raw_brenda = Path(config.project_root) / "data" / "raw" / "brenda"
    raw_metacyc = Path(config.project_root) / "data" / "raw" / "metacyc"
    raw_mapping_ec_numbers = Path(config.project_root) / "data" / "raw" / "mapping_ec_numbers"
    all_ok &= check_path_exists(raw_brenda, "BRenda raw directory")
    all_ok &= check_path_exists(raw_metacyc, "MetaCyc raw directory")
    all_ok &= check_path_exists(raw_mapping_ec_numbers, "mapping_ec_numbers raw directory")

    # Check the dynamically discovered BRenda JSON file.
    if hasattr(config, "brenda_json_path"):
        brenda_path = Path(config.brenda_json_path)
        if brenda_path.exists() and brenda_path.is_file() and brenda_path.suffix.lower() == ".json":
            print(f"[OK] BRenda JSON file exists and is valid: {brenda_path}")
        else:
            print(f"[ERROR] BRenda JSON file is invalid, not a file, or does not have a .json extension: {brenda_path}")
            all_ok = False
    else:
        print("[ERROR] Config does not define 'brenda_json_path'.")
        all_ok = False

    # Check the dynamically discovered MetaCyc directory.
    if hasattr(config, "metacyc_directory_path"):
        metacyc_path = Path(config.metacyc_directory_path)
        if not metacyc_path.exists():
            print(f"[ERROR] MetaCyc path does not exist: {metacyc_path}")
            all_ok = False
        elif metacyc_path.is_file():
            print(f"[ERROR] MetaCyc path is a file, but should be a directory: {metacyc_path}")
            all_ok = False
        elif metacyc_path.is_dir():
            # Check if the directory name contains a dot.
            if '.' in metacyc_path.name:
                print(f"[ERROR] MetaCyc directory name should not contain a period: {metacyc_path.name}")
                all_ok = False
            else:
                print(f"[OK] MetaCyc directory exists: {metacyc_path}")
        else:
            print(f"[ERROR] MetaCyc path is neither a file nor a directory: {metacyc_path}")
            all_ok = False
    else:
        print("[ERROR] Config does not define 'metacyc_directory_path'.")
        all_ok = False

    # For output file paths, check that their parent directories exist.
    output_checks = [
        ("brenda_extracted_ids_output_path", getattr(config, "brenda_extracted_ids_output_path", None)),
        ("brenda_reaction_csv", getattr(config, "brenda_reaction_csv", None)),
        ("brenda_kcat_csv", getattr(config, "brenda_kcat_csv", None)),
        ("kegg_extracted_ids_output_path", getattr(config, "kegg_extracted_ids_output_path", None)),
        ("kegg_reaction_csv", getattr(config, "kegg_reaction_csv", None)),
        ("kegg_reaction_pathway_csv", getattr(config, "kegg_reaction_pathway_csv", None)),
        ("metacyc_extracted_ids_output_path", getattr(config, "metacyc_extracted_ids_output_path", None)),
        ("metacyc_reaction_csv", getattr(config, "metacyc_reaction_csv", None)),
        ("metacyc_kcat_csv", getattr(config, "metacyc_kcat_csv", None)),
        ("metacyc_reaction_pathway_csv", getattr(config, "metacyc_reaction_pathway_csv", None)),
        ("metacyc_pathway_csv", getattr(config, "metacyc_pathway_csv", None)),
        ("metacyc_reaction_csv_with_smiles", getattr(config, "metacyc_reaction_csv_with_smiles", None)),
        ("merged_extracted_ids_output_path", getattr(config, "merged_extracted_ids_output_path", None)),
        ("merged_ec_numbers_with_kcat", getattr(config, "merged_ec_numbers_with_kcat", None)),
        ("sabio_temp_data_ec_number", getattr(config, "sabio_temp_data_ec_number", None)),
        ("sabio_temp_data_kcat_unisubstrate", getattr(config, "sabio_temp_data_kcat_unisubstrate", None)),
        ("sabio_kcat_csv", getattr(config, "sabio_kcat_csv", None)),
        ("predictkcat_input_tsv", getattr(config, "predictkcat_input_tsv", None)),
        ("predictkcat_output_tsv", getattr(config, "predictkcat_output_tsv", None)),
        ("predictkcat_input_with_ec", getattr(config, "predictkcat_input_with_ec", None)),
        ("predictkcat_with_ec_protein_sequences", getattr(config, "predictkcat_with_ec_protein_sequences", None)),
        ("predictkcat_missing_kcat_ids_output_path", getattr(config, "predictkcat_missing_kcat_ids_output_path", None)),
        ("predictkcat_failed_kcat_ids_output_path", getattr(config, "predictkcat_failed_kcat_ids_output_path", None)),
        ("ecs_with_kcat_output_path", getattr(config, "ecs_with_kcat_output_path", None)),
        ("mapping_ec_number_old_new", getattr(config, "mapping_ec_number_old_new", None)),
        ("merged_overview", getattr(config, "merged_overview", None)),
    ]
    for desc, path in output_checks:
        if path:
            all_ok &= check_parent_exists(path, desc)
        else:
            print(f"[WARNING] Config does not define '{desc}'.")

    if all_ok:
        print("\nAll prechecks passed successfully.")
        sys.exit(0)
    else:
        print("\nPrecheck failed. Please fix the issues above and try again.")
        sys.exit(1)

if __name__ == "__main__":
    main()
