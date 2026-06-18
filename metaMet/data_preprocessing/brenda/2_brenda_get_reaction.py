import pandas as pd
import ijson

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

# Read EC numbers from a text file into a list
with open(config.brenda_extracted_ids_output_path, 'r') as ec_file:
    ec_numbers = [line.strip() for line in ec_file if line.strip()]

records = []

# Loop through the EC numbers
for ec_number in ec_numbers:
    print(ec_number)
    
    # Using a set to track unique reactions
    reaction_set = set()
    
    # Open the JSON file again for each EC number to reset the file pointer
    with open(config.brenda_json_path, 'r') as json_file:
        ec_path = f'data.{ec_number}'
        try:
            # Use ijson to parse the file incrementally for the given EC number
            for item in ijson.items(json_file, ec_path):
                if 'reaction' in item:
                    for reaction in item['reaction']:
                        educts = tuple(reaction.get('educts', []))
                        products = tuple(reaction.get('products', []))
                        
                        # Only add reactions without '?' and ensure uniqueness by using a set
                        if '?' not in educts and '?' not in products:
                            reaction_set.add((ec_number, 'reaction', educts, products))
        
        except ijson.JSONError as e:
            print(f"Error processing EC number {ec_number}: {e}")
        except StopIteration:
            print(f"No data found for EC number {ec_number}")
    
    # Convert set items back to records for further processing
    for ec, type_, educts, products in reaction_set:
        records.append({'ec_number': ec, 'type': type_, 'educts': list(educts), 'products': list(products)})

# Convert the collected records to a DataFrame
df = pd.DataFrame(records)

# Write the DataFrame to a CSV file
df.to_csv(config.brenda_reaction_csv, index=False)

print("Reactions data has been successfully written to 'brenda_reactions.csv'.")

# Check for 'generic_reaction', 'natural_reaction', or 'reaction'
# if 'generic_reaction' in item:
#     for reaction in item['generic_reaction']:
#         educts = reaction.get('educts', [])
#         products = reaction.get('products', [])
#         records.append({'ec_number': ec_number, 'educts': educts, 'products': products})

# If you want to uncomment these parts, it will handle 'natural_reaction' and 'reaction'
# if 'natural_reaction' in item:
#     for reaction in item['natural_reaction']:
#         educts = reaction.get('educts', [])
#         products = reaction.get('products', [])
#         records.append({'ec_number': ec_number, 'type': 'natural_reaction', 'educts': educts, 'products': products})


