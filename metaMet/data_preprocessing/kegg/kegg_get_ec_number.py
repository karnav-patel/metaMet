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

def fetch_ec_numbers_from_kegg():
    # Fetch all EC numbers available in KEGG
    url = "http://rest.kegg.jp/list/enzyme"
    response = requests.get(url)
    
    if response.status_code == 200:
        # Extract EC numbers from the response
        ec_numbers = [line.split("\t")[0].replace("ec:", "") for line in response.text.splitlines()]
        return ec_numbers
    else:
        print(f"Error fetching EC numbers: {response.status_code}")
        return []

def save_ec_numbers_to_file(ec_numbers, filename):
    # Save EC numbers to a text file
    with open(filename, 'w') as file:
        for ec_number in ec_numbers:
            file.write(f"{ec_number}\n")

# Fetch EC numbers from KEGG
ec_numbers = fetch_ec_numbers_from_kegg()

# Save EC numbers to a file
save_ec_numbers_to_file(ec_numbers, config.kegg_extracted_ids_output_path)

print(f"Saved {len(ec_numbers)} EC numbers to ec_numbers.txt")
