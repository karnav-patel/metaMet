import json

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


def extract_enzyme_codes_recursive(data, enzyme_codes=None):
    """
    Recursively search through nested dictionaries and lists to extract enzyme codes.
    An enzyme code is identified as a key with more than one dot where each part is numeric.
    """
    if enzyme_codes is None:
        enzyme_codes = set()  # Use a set to prevent duplicates

    if isinstance(data, dict):
        for key, value in data.items():
            # Check if the key resembles an enzyme code (e.g., "1.1.1")
            if key.count('.') > 1:
                parts = key.split('.')
                if all(part.isdigit() for part in parts):
                    enzyme_codes.add(key)
            # Recursively process the value if it's a dict or list
            extract_enzyme_codes_recursive(value, enzyme_codes)
    elif isinstance(data, list):
        for item in data:
            extract_enzyme_codes_recursive(item, enzyme_codes)

    return enzyme_codes


def main():
    # Load JSON data using the file path from config
    with open(config.brenda_json_path, 'r', encoding='utf-8') as f:
        brenda_data = json.load(f)

    # Extract enzyme codes recursively
    enzyme_codes = extract_enzyme_codes_recursive(brenda_data)

    # Save the extracted enzyme codes to the output file path from config
    with open(config.brenda_extracted_ids_output_path, 'w', encoding='utf-8') as output_file:
        for code in sorted(enzyme_codes):
            output_file.write(code + '\n')

    print(f'Extracted {len(enzyme_codes)} unique enzyme codes.')
    print(f'Enzyme codes have been saved to {config.brenda_extracted_ids_output_path}')


if __name__ == '__main__':
    main()
