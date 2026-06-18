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
    """
    Fetch reaction IDs associated with the given EC number from KEGG.
    """
    url = f"http://rest.kegg.jp/link/reaction/enzyme:{ec_number}"
    response = requests.get(url)

    if response.status_code == 200:
        reactions = []
        for line in response.text.splitlines():
            parts = line.split("\t")
            if len(parts) >= 2:
                reactions.append(f"'{parts[1].strip()}'")  # Store reactions in proper list format
        return reactions
    else:
        print(f"Error fetching reactions for EC {ec_number}: {response.status_code}")
        return []

def fetch_pathways_by_reaction(reaction_ids):
    """
    Fetch pathway IDs associated with the given list of reaction IDs.
    """
    pathway_set = set()

    for reaction_id in reaction_ids:
        reaction_id = reaction_id.strip("'")  # Remove single quotes for URL construction
        url = f"http://rest.kegg.jp/link/pathway/{reaction_id}"
        response = requests.get(url)

        if response.status_code == 200:
            for line in response.text.splitlines():
                parts = line.split("\t")
                if len(parts) >= 2:
                    pathway_set.add(f"'{parts[1].strip()}'")  # Store pathways in proper list format
        else:
            print(f"Error fetching pathways for reaction {reaction_id}: {response.status_code}")

    return list(pathway_set)

def process_ec_numbers(input_file, output_file):
    """
    Read EC numbers from the input file, retrieve reaction IDs and pathway IDs,
    and write the results to a CSV file.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["ec_number", "reaction_ids", "pathway_ids"])

        for line in infile:
            ec_number = line.strip()
            if ec_number:
                print(f"Processing EC: {ec_number}")

                # Fetch reactions and pathways
                reaction_ids = fetch_reaction_by_ec(ec_number)
                pathway_ids = fetch_pathways_by_reaction(reaction_ids)

                # Convert lists to properly formatted strings
                reaction_str = f"[{', '.join(reaction_ids)}]"
                pathway_str = f"[{', '.join(pathway_ids)}]"

                # Write to CSV
                csv_writer.writerow([ec_number, reaction_str, pathway_str])

# Usage example
# Usage example
input_file = config.kegg_extracted_ids_output_path  # File containing EC numbers, one per line
output_file = config.kegg_reaction_pathway_csv  # Output CSV file

process_ec_numbers(input_file, output_file)
