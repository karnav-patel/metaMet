import csv
import requests

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

def fetch_reaction_by_ec(ec_number):
    # Fetch reaction IDs associated with the EC number
    url = f"http://rest.kegg.jp/link/reaction/enzyme:{ec_number}"
    response = requests.get(url)

    if response.status_code == 200:
        reactions = []
        for line in response.text.splitlines():
            parts = line.split("\t")
            if len(parts) >= 2:
                reactions.append(parts[1].strip())
        return reactions
    else:
        print(f"Error fetching reactions for EC {ec_number}: {response.status_code}")
        return []

def fetch_reaction_details(reaction_id):
    # Fetch reaction details by reaction ID
    url = f"http://rest.kegg.jp/get/{reaction_id}"
    response = requests.get(url)

    if response.status_code == 200:
        return response.text
    else:
        print(f"Error fetching details for reaction {reaction_id}: {response.status_code}")
        return ""

def extract_reactants_products(details_text):
    # Extract reactants and products from the definition line
    lines = details_text.splitlines()
    for line in lines:
        if line.startswith("DEFINITION"):
            definition = line.split("DEFINITION")[1].strip()
            if "<=>" in definition:
                reactants, products = definition.split("<=>")
            elif "=" in definition:
                reactants, products = definition.split("=")
            else:
                return [], []

            # Split reactants and products into lists
            reactants = [r.strip() for r in reactants.split(" + ")]
            products = [p.strip() for p in products.split(" + ")]
            return reactants, products
    return [], []

def process_ec_numbers(input_file, output_file):
    # Read EC numbers from the input file and write results to a CSV file
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["ec_number", "educts", "products"])

        for line in infile:
            ec_number = line.strip()
            if ec_number:
                print(f"Processing EC: {ec_number}")
                reactions = fetch_reaction_by_ec(ec_number)
                for reaction_id in reactions:
                    details = fetch_reaction_details(reaction_id)
                    reactants, products = extract_reactants_products(details)
                    if reactants and products:
                        csv_writer.writerow([ec_number, reactants, products])

# Usage example
input_file = config.kegg_extracted_ids_output_path  # File containing EC numbers, one per line
output_file = config.kegg_reaction_csv  # Output CSV file
process_ec_numbers(input_file, output_file)
