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

def clean_html(name):
    """
    Cleans HTML from the given name by:
    - Replacing <sup>+</sup> with +
    - Removing <i> tags while keeping the text inside
    - Decoding HTML entities
    - Removing any other HTML tags
    """
    if not name:
        return name

    # Replace <sup>+</sup> with regular plus
    name = re.sub(r'<sup>\+</sup>', '+', name, flags=re.IGNORECASE)

    # Replace <i>text</i> with text (remove italics)
    name = re.sub(r'<i>(.*?)</i>', r'\1', name, flags=re.IGNORECASE)

    # Replace common HTML entities with Unicode characters using html.unescape
    name = html.unescape(name)

    # Remove any other remaining HTML tags
    name = re.sub(r'<[^>]+>', '', name)

    return name

def parse_pathway_data(pathway_file):
    """
    Parses pathway information from the pathways.dat file.

    Parameters:
    - pathway_file: Path to the pathways.dat file.

    Returns:
    - A list of dictionaries containing PathwayID, Common-Name, Species, Sub-Pathways, and Super-Pathways.
    """
    pathway_entries = []
    pathway_id = None
    common_name = None
    species_list = []
    sub_pathway_list = []
    super_pathway_list = []

    with codecs.open(pathway_file, 'r', encoding='utf-8', errors='ignore') as file:
        for line in file:
            line = line.strip()

            # Detect the end of a pathway entry
            if line == "//":
                if pathway_id:
                    pathway_entries.append({
                        "PathwayID": pathway_id,
                        "Common-Name": clean_html(common_name),  # Clean HTML before storing
                        "Species": str(species_list),
                        "Sub-Pathways": str(sub_pathway_list),
                        "Super-Pathways": str(super_pathway_list)
                    })
                # Reset for the next entry
                pathway_id = None
                common_name = None
                species_list = []
                sub_pathway_list = []
                super_pathway_list = []
                continue

            # Parse fields
            if line.startswith("UNIQUE-ID"):
                pathway_id = line.split(" - ", 1)[1].strip()
            elif line.startswith("COMMON-NAME"):
                common_name = line.split(" - ", 1)[1].strip()
            elif line.startswith("SPECIES"):
                species_list.append(line.split(" - ", 1)[1].strip())
            elif line.startswith("SUB-PATHWAYS"):
                sub_pathway_list.append(line.split(" - ", 1)[1].strip())
            elif line.startswith("SUPER-PATHWAYS"):
                super_pathway_list.append(line.split(" - ", 1)[1].strip())

    return pathway_entries

def search_and_extract_pathways(root_dir, output_csv):
    """
    Searches for pathways.dat files in the specified root directory,
    extracts PathwayID, Common-Name, Species, Sub-Pathways, and Super-Pathways,
    and writes them to a CSV file.

    Parameters:
    - root_dir: Root directory to start the search.
    - output_csv: Path to the output CSV file where pathway data will be saved.
    """
    all_pathways = []

    # Walk through each directory and file in the root_dir to parse all pathways.dat files
    for dirpath, _, filenames in os.walk(root_dir):
        if 'pathways.dat' in filenames:
            pathway_file_path = os.path.join(dirpath, 'pathways.dat')
            print(f"Processing pathway file: {pathway_file_path}")
            pathways = parse_pathway_data(pathway_file_path)
            all_pathways.extend(pathways)

    # Write data to CSV
    try:
        with open(output_csv, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['PathwayID', 'Common-Name', 'Species', 'Sub-Pathways', 'Super-Pathways']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for pathway in all_pathways:
                writer.writerow(pathway)
        print(f"Successfully saved pathway data to '{output_csv}'.")
    except Exception as e:
        print(f"An unexpected error occurred while writing to the CSV file: {e}")

# Specify the root directory and output CSV file path
root_directory = config.metacyc_directory_path
output_csv = config.metacyc_pathway_csv

# Call the function to search for pathways.dat files,
# extract PathwayID, Common-Name, Species, Sub-Pathways, Super-Pathways, and save to CSV
search_and_extract_pathways(root_directory, output_csv)
