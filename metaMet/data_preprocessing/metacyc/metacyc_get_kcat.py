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

def clean_compound_name(name):
    """
    Removes HTML tags from the compound name and replaces HTML entities with Unicode characters.
    """
    if not name:
        return name

    # Replace <sup>+</sup> with plus
    name = re.sub(r'<sup>\+</sup>', '+', name, flags=re.IGNORECASE)
    # Remove <i>...</i> tags
    name = re.sub(r'<i>(.*?)</i>', r'\1', name, flags=re.IGNORECASE)
    # Decode HTML entities
    name = html.unescape(name)
    # Remove any remaining HTML tags
    clean = re.sub(r'<[^>]+>', '', name)
    return clean

def parse_compounds_data(compounds_file, compounds_mapping):
    unique_id = None
    common_name = None

    with codecs.open(compounds_file, 'r', encoding='utf-8', errors='ignore') as file:
        for line in file:
            line = line.strip()
            if line == "//":
                if unique_id and common_name:
                    compounds_mapping[unique_id] = common_name
                unique_id = None
                common_name = None
                continue

            if line.startswith("UNIQUE-ID"):
                unique_id = line.split(" - ", 1)[1].strip()
            elif line.startswith("COMMON-NAME"):
                common_name = line.split(" - ", 1)[1].strip()

def parse_enzrxns_data(enzrxns_file, reaction_kcat_map):
    block_reaction_ids = []
    block_pairs = []
    current_kcat = None

    with codecs.open(enzrxns_file, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if line == "//":
                for rxn_id in block_reaction_ids:
                    reaction_kcat_map[rxn_id].extend(block_pairs)
                block_reaction_ids = []
                block_pairs = []
                current_kcat = None
                continue

            if not line:
                continue

            if line.startswith("UNIQUE-ID -"):
                continue

            if line.startswith("REACTION -"):
                reaction_id = line.split(" - ", 1)[1].strip()
                block_reaction_ids.append(reaction_id)
            elif line.startswith("KCAT -"):
                kcat_str = line.split(" - ", 1)[1].strip()
                try:
                    current_kcat = float(kcat_str)
                except ValueError:
                    current_kcat = None
            elif line.startswith("^SUBSTRATE -"):
                if current_kcat is not None:
                    substrate_id = line.split(" - ", 1)[1].strip()
                    block_pairs.append((current_kcat, substrate_id))
                    current_kcat = None
            else:
                if not line.startswith("^SUBSTRATE"):
                    current_kcat = None

def parse_reaction_data(input_file, reaction_data, compounds_mapping, reaction_kcat_map):
    unique_id = None
    ec_numbers = []
    reactants = []
    products = []

    with codecs.open(input_file, 'r', encoding='utf-8', errors='ignore') as file:
        for line in file:
            line = line.strip()
            if line == "//":
                if unique_id:
                    reactant_names = [
                        clean_compound_name(compounds_mapping.get(cpd, cpd))
                        for cpd in reactants
                    ]
                    product_names = [
                        clean_compound_name(compounds_mapping.get(cpd, cpd))
                        for cpd in products
                    ]

                    pairs = reaction_kcat_map.get(unique_id, [])
                    if pairs:
                        for kcat, substrate_id in pairs:
                            substrate_name = clean_compound_name(
                                compounds_mapping.get(substrate_id, substrate_id)
                            )
                            reaction_data.append({
                                "ReactionID": unique_id,
                                "ec_number": ", ".join(ec_numbers),
                                "educts": str(reactant_names),
                                "products": str(product_names),
                                "kcat": kcat,
                                "substrate": substrate_name
                            })
                    else:
                        reaction_data.append({
                            "ReactionID": unique_id,
                            "ec_number": ", ".join(ec_numbers),
                            "educts": str(reactant_names),
                            "products": str(product_names),
                            "kcat": "",
                            "substrate": ""
                        })
                unique_id = None
                ec_numbers = []
                reactants = []
                products = []
                continue

            if line.startswith("UNIQUE-ID"):
                unique_id = line.split(" - ", 1)[1].strip()
            elif line.startswith("EC-NUMBER"):
                ec_number = line.split(" - ", 1)[1].strip().lstrip("EC-")
                ec_numbers.append(ec_number)
            elif line.startswith("LEFT"):
                reactants.append(line.split(" - ", 1)[1].strip())
            elif line.startswith("RIGHT"):
                products.append(line.split(" - ", 1)[1].strip())

def search_and_extract_reaction_data(root_dir, output_csv):
    reaction_data = []
    compounds_mapping = {}
    reaction_kcat_map = defaultdict(list)

    # 1) compounds.dat
    for dirpath, _, filenames in os.walk(root_dir):
        if 'compounds.dat' in filenames:
            parse_compounds_data(os.path.join(dirpath, 'compounds.dat'), compounds_mapping)

    # 2) enzrxns.dat
    for dirpath, _, filenames in os.walk(root_dir):
        if 'enzrxns.dat' in filenames:
            parse_enzrxns_data(os.path.join(dirpath, 'enzrxns.dat'), reaction_kcat_map)

    # 3) reactions.dat
    for dirpath, _, filenames in os.walk(root_dir):
        if 'reactions.dat' in filenames:
            parse_reaction_data(os.path.join(dirpath, 'reactions.dat'),
                                reaction_data,
                                compounds_mapping,
                                reaction_kcat_map)

    # de-dup & sort
    unique_rows = []
    seen = set()
    for row in reaction_data:
        key = (row["ReactionID"], row["ec_number"], row["educts"],
               row["products"], row["kcat"], row["substrate"])
        if key not in seen:
            seen.add(key)
            unique_rows.append(row)

    def sorting_key(r):
        if r["kcat"] != "":
            try:
                kv = float(r["kcat"])
            except:
                kv = 0
            return (0, -kv)
        return (1, 0)

    unique_rows.sort(key=sorting_key)

    # write CSV
    try:
        with open(output_csv, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['ReactionID', 'ec_number', 'educts', 'products', 'kcat', 'substrate']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for data in unique_rows:
                writer.writerow(data)
        print(f"Saved reaction data to '{output_csv}'.")
    except Exception as e:
        print(f"Error writing CSV: {e}")

if __name__ == "__main__":
    root_directory = config.metacyc_directory_path
    output_csv      = config.metacyc_kcat_csv

    search_and_extract_reaction_data(root_directory, output_csv)
