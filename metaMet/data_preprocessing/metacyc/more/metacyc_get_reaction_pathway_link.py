import os
import csv
import codecs

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

def parse_reactions_data(reactions_file):
    """
    Parses reaction information from the reactions.dat file.

    Parameters:
    - reactions_file: Path to the reactions.dat file.

    Returns:
    - A list of dictionaries containing UNIQUE-ID, PathwayIDs, and OtherReactionIds.
    """
    reaction_entries = []
    unique_id = None
    pathway_ids = []
    other_reaction_ids = []

    with codecs.open(reactions_file, 'r', encoding='utf-8', errors='ignore') as file:
        for line in file:
            line = line.strip()

            # Detect the end of a reaction entry
            if line == "//":
                if unique_id:
                    reaction_entries.append({
                        "UNIQUE-ID": unique_id,
                        "PathwayIDs": pathway_ids.copy(),
                        "OtherReactionIds": other_reaction_ids.copy()
                    })
                # Reset for the next reaction
                unique_id = None
                pathway_ids = []
                other_reaction_ids = []
                continue

            # Parse fields
            if line.startswith("UNIQUE-ID"):
                unique_id = line.split(" - ", 1)[1].strip()
            elif line.startswith("IN-PATHWAY"):
                pathway_id = line.split(" - ", 1)[1].strip()
                pathway_ids.append(pathway_id)
            elif line.startswith("REACTION-LIST"):
                reaction_id = line.split(" - ", 1)[1].strip()
                other_reaction_ids.append(reaction_id)

    return reaction_entries

def search_and_extract_reactions(root_dir, output_csv):
    """
    Searches for reactions.dat files in the specified root directory,
    extracts UNIQUE-ID, PathwayIDs, OtherReactionIds, and the directory name,
    removes duplicate rows (if the entire row is a duplicate),
    and writes the unique data to a CSV file.

    Parameters:
    - root_dir: Root directory to start the search.
    - output_csv: Path to the output CSV file where reaction data will be saved.
    """
    all_reactions = []

    # Walk through each directory and file in the root_dir to parse all reactions.dat files
    for dirpath, _, filenames in os.walk(root_dir):
        if 'reactions.dat' in filenames:
            reactions_file_path = os.path.join(dirpath, 'reactions.dat')
            print(f"Processing reactions file: {reactions_file_path}")
            reactions = parse_reactions_data(reactions_file_path)
            # Get the directory name (last part of the dirpath)
            directory_name = os.path.basename(dirpath)
            for reaction in reactions:
                reaction['DirectoryName'] = directory_name
            all_reactions.extend(reactions)

    # Remove duplicate rows (if the whole row matches)
    unique_reactions = []
    seen = set()
    for reaction in all_reactions:
        # Convert lists to tuples so they are hashable for the set, and include DirectoryName
        reaction_tuple = (
            reaction['UNIQUE-ID'],
            tuple(reaction['PathwayIDs']),
            tuple(reaction['OtherReactionIds']),
            reaction['DirectoryName']
        )
        if reaction_tuple not in seen:
            seen.add(reaction_tuple)
            unique_reactions.append(reaction)

    # Write unique data to CSV
    try:
        with open(output_csv, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['UNIQUE-ID', 'PathwayIDs', 'OtherReactionIds', 'DirectoryName']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for reaction in unique_reactions:
                writer.writerow({
                    'UNIQUE-ID': reaction['UNIQUE-ID'],
                    'PathwayIDs': str(reaction['PathwayIDs']),
                    'OtherReactionIds': str(reaction['OtherReactionIds']),
                    'DirectoryName': reaction['DirectoryName']
                })
        print(f"Successfully saved reaction data to '{output_csv}'.")
    except Exception as e:
        print(f"An unexpected error occurred while writing to the CSV file: {e}")

# Specify the root directory and output CSV file path
root_directory = config.metacyc_directory_path
output_csv = config.metacyc_reaction_pathway_csv

# Call the function to search for reactions.dat files,
# extract UNIQUE-ID, PathwayIDs, OtherReactionIds, add DirectoryName,
# remove duplicates, and save to CSV
search_and_extract_reactions(root_directory, output_csv)
