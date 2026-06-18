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

def extract_unique_ec_numbers(input_file, ec_numbers_set):
    """
    Extracts unique EC numbers from the input file and adds them to a set.
    
    Parameters:
    - input_file: Path to the reactions.dat file.
    - ec_numbers_set: Set to store unique EC numbers.
    """
    ec_number_pattern = re.compile(r'^EC-NUMBER\s*-\s*EC-([\d\.]+)', re.IGNORECASE)

    try:
        # Attempt to open the file with UTF-8 encoding first
        with open(input_file, 'r', encoding='utf-8') as file:
            lines = file.readlines()
    except UnicodeDecodeError:
        # Fallback to ISO-8859-1 if UTF-8 fails
        try:
            print(f"Retrying with ISO-8859-1 encoding for file: {input_file}")
            with open(input_file, 'r', encoding='ISO-8859-1') as file:
                lines = file.readlines()
        except Exception as fallback_e:
            print(f"Error reading file '{input_file}' with fallback encoding: {fallback_e}")
            return
    except FileNotFoundError:
        print(f"Error: The file '{input_file}' does not exist.")
        return
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return

    # Process each line for EC numbers
    for line in lines:
        ec_match = ec_number_pattern.match(line.strip())
        if ec_match:
            ec_numbers_set.add(ec_match.group(1))

def search_and_extract_ec_numbers(root_dir, output_file):
    """
    Searches for reactions.dat files in the specified root directory, extracts
    unique EC numbers, sorts them, and writes them to the output file.
    
    Parameters:
    - root_dir: Root directory to start the search.
    - output_file: Path to the output text file where sorted EC numbers will be saved.
    """
    ec_numbers_set = set()

    # Walk through each directory and file in the root_dir
    for dirpath, _, filenames in os.walk(root_dir):
        if 'reactions.dat' in filenames:
            input_file_path = os.path.join(dirpath, 'reactions.dat')
            print(f"Processing file: {input_file_path}")
            extract_unique_ec_numbers(input_file_path, ec_numbers_set)

    # Sort EC numbers and write to output file
    sorted_ec_numbers = sorted(ec_numbers_set)

    try:
        with open(output_file, 'w', encoding='utf-8') as txtfile:
            for ec_number in sorted_ec_numbers:
                txtfile.write(ec_number + "\n")
        print(f"Successfully saved sorted EC numbers to '{output_file}'.")
    except Exception as e:
        print(f"An unexpected error occurred while writing to the file: {e}")

# Specify the root directory and output file path
root_directory = config.metacyc_directory_path
output_filepath = config.metacyc_extracted_ids_output_path

# Call the function to search for reactions.dat files and extract EC numbers
search_and_extract_ec_numbers(root_directory, output_filepath)
