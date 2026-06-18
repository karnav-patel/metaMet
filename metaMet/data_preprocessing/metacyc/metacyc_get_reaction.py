import os
import csv
import codecs
import re
import html
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

def clean_compound_name(name):
    if not name:
        return name
    name = re.sub(r'<sup>\+</sup>', '+', name, flags=re.IGNORECASE)
    name = re.sub(r'<i>(.*?)</i>', r'\1', name, flags=re.IGNORECASE)
    name = html.unescape(name)
    return re.sub(r'<[^>]+>', '', name)

def parse_compounds_data(compounds_file, name_map, smiles_map):
    unique_id = None
    with codecs.open(compounds_file, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if line == "//":
                unique_id = None
                continue
            if line.startswith("UNIQUE-ID"):
                unique_id = line.split(" - ",1)[1].strip()
            elif line.startswith("COMMON-NAME") and unique_id:
                name_map[unique_id] = line.split(" - ",1)[1].strip()
            elif line.startswith("SMILES") and unique_id:
                smiles_map[unique_id] = line.split(" - ",1)[1].strip()

def parse_reaction_data(input_file, reaction_data, name_map, smiles_map):
    unique_id = None
    ec_numbers = []
    reactants = []
    products = []

    with codecs.open(input_file, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if line == "//":
                if unique_id:
                    # build cleaned names
                    reactant_names = [ clean_compound_name(name_map.get(r,r)) for r in reactants ]
                    product_names = [ clean_compound_name(name_map.get(p,p)) for p in products ]
                    # build SMILES lists
                    reactant_smiles = [ smiles_map.get(r, "") for r in reactants ]
                    product_smiles = [ smiles_map.get(p, "") for p in products ]

                    reaction_data.append({
                        "ReactionID": unique_id,
                        "ec_number": ", ".join(ec_numbers),
                        "educts": str(reactant_names),
                        "products": str(product_names),
                        "educts_smiles": str(reactant_smiles),
                        "products_smiles": str(product_smiles),
                    })
                # reset
                unique_id = None
                ec_numbers = []
                reactants = []
                products = []
                continue

            if line.startswith("UNIQUE-ID"):
                unique_id = line.split(" - ",1)[1].strip()
            elif line.startswith("EC-NUMBER"):
                ec = line.split(" - ",1)[1].strip().lstrip("EC-")
                ec_numbers.append(ec)
            elif line.startswith("LEFT") and unique_id:
                reactants.append(line.split(" - ",1)[1].strip())
            elif line.startswith("RIGHT") and unique_id:
                products.append(line.split(" - ",1)[1].strip())

def search_and_extract_reaction_data(root_dir, output_csv, output_csv_with_smiles):
    reaction_data = []
    name_map = {}
    smiles_map = {}

    # parse all compounds first
    for dirpath, _, files in os.walk(root_dir):
        if 'compounds.dat' in files:
            parse_compounds_data(os.path.join(dirpath,'compounds.dat'),
                                 name_map, smiles_map)

    # parse all reactions
    for dirpath, _, files in os.walk(root_dir):
        if 'reactions.dat' in files:
            parse_reaction_data(os.path.join(dirpath,'reactions.dat'),
                                reaction_data, name_map, smiles_map)

    # sorting helpers (unchanged)
    def ec_sort_key(ec):
        if ec == "1.1.1.1": return (0,)
        if ec:
            try: return (1, tuple(int(p) for p in ec.split('.')))
            except: return (1,())
        return (2,)
    def reaction_sort_key(r):
        first = r['ec_number'].split(", ")[0] if r['ec_number'] else ''
        return ec_sort_key(first)

    reaction_data.sort(key=reaction_sort_key)

    # write original CSV
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['ReactionID','ec_number','educts','products'])
        writer.writeheader()
        for r in reaction_data:
            writer.writerow({k:r[k] for k in writer.fieldnames})

    # write CSV *with* SMILES
    with open(output_csv_with_smiles, 'w', newline='', encoding='utf-8') as f:
        fn = ['ReactionID','ec_number','educts','products','educts_smiles','products_smiles']
        writer = csv.DictWriter(f, fieldnames=fn)
        writer.writeheader()
        for r in reaction_data:
            writer.writerow({k:r[k] for k in fn})

if __name__ == "__main__":
    root_directory            = config.metacyc_directory_path
    output_csv                = config.metacyc_reaction_csv
    output_csv_with_smiles    = config.metacyc_reaction_csv_with_smiles

    search_and_extract_reaction_data(
        root_directory,
        output_csv,
        output_csv_with_smiles
    )

