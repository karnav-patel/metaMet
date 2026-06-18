from pathlib import Path

# Get the current file's absolute path.
current_file_path = Path(__file__).resolve()

# Traverse upwards until we find the "metaMet" folder.
for parent in current_file_path.parents:
    if parent.name == "metaMet":
        project_root = parent
        break
else:
    raise RuntimeError("Project root folder 'metaMet' not found.")

# --- Dynamically locate the single BRenda JSON file ---
brenda_dir = project_root / "data" / "raw" / "brenda"
# List only files (ignore directories)
brenda_files = [f for f in brenda_dir.iterdir() if f.is_file()]
if len(brenda_files) != 1:
    raise RuntimeError(
        f"Expected exactly one file in {brenda_dir}, but found {len(brenda_files)}."
    )
brenda_json_path = brenda_files[0]

# --- Dynamically locate the single MetaCyc folder ---
metacyc_raw_dir = project_root / "data" / "raw" / "metacyc"
# List only directories
metacyc_subdirs = [d for d in metacyc_raw_dir.iterdir() if d.is_dir()]
if len(metacyc_subdirs) != 1:
    raise RuntimeError(
        f"Expected exactly one folder in {metacyc_raw_dir}, but found {len(metacyc_subdirs)}."
    )
metacyc_directory_path = metacyc_subdirs[0]

# --- Define the rest of your file paths relative to project_root ---
brenda_extracted_ids_output_path = project_root / "data" / "processed" / "brenda" / "brenda_processed_ec_numbers.txt"
brenda_reaction_csv = project_root / "data" / "processed" / "brenda" / "brenda_reactions.csv"
brenda_kcat_csv = project_root / "data" / "processed" / "brenda" / "brenda_kcat.csv"

kegg_extracted_ids_output_path = project_root / "data" / "processed" / "kegg" / "kegg_processed_ec_numbers.txt"
kegg_reaction_csv = project_root / "data" / "processed" / "kegg" / "kegg_reactions.csv"
kegg_reaction_pathway_csv = project_root / "data" / "processed" / "kegg" / "more" / "kegg_reaction_pathway_link.csv"

metacyc_extracted_ids_output_path = project_root / "data" / "processed" / "metacyc" / "metacyc_processed_ec_numbers.txt"
metacyc_reaction_csv = project_root / "data" / "processed" / "metacyc" / "metacyc_reactions.csv"
metacyc_kcat_csv = project_root / "data" / "processed" / "metacyc" / "metacyc_kcat.csv"
metacyc_reaction_pathway_csv = project_root / "data" / "processed" / "metacyc" / "more" / "metacyc_reaction_pathway_link.csv"
metacyc_pathway_csv = project_root / "data" / "processed" / "metacyc" / "more" / "metacyc_pathway.csv"
metacyc_reaction_csv_with_smiles = project_root / "data" / "processed" / "metacyc" / "more" / "metacyc_reactions_with_smiles.csv"

merged_extracted_ids_output_path = project_root / "data" / "processed" / "merged" / "all_processed_ec_numbers.txt"
merged_ec_numbers_with_kcat = project_root / "data" / "processed" / "merged" / "all_processed_ec_numbers_with_kcat.txt"

sabio_temp_data_ec_number = project_root / "data" / "processed" / "sabio" / "more" / "ec"
sabio_temp_data_kcat_unisubstrate = project_root / "data" / "processed" / "sabio" / "more" / "sabio_kcat_unisubstrate.tsv"
sabio_kcat_csv = project_root / "data" / "processed" / "sabio" / "sabio_kcat.tsv"

predictkcat_input_tsv = project_root / "data" / "processed" / "predictkcat" / "input.tsv"
predictkcat_output_tsv = project_root / "data" / "processed" / "predictkcat" / "output.tsv"
predictkcat_input_with_ec = project_root / "data" / "processed" / "predictkcat" / "input_with_ec.csv"
predictkcat_with_ec_protein_sequences = project_root / "data" / "processed" / "predictkcat" / "predictkcat_with_ec_protein_sequences.csv"
predictkcat_missing_kcat_ids_output_path = project_root / "data" / "processed" / "predictkcat" / "predictkcat_missing_kcat_ids_output_path.txt"
predictkcat_failed_kcat_ids_output_path = project_root / "data" / "processed" / "predictkcat" / "predictkcat_failed_kcat_ids_output_path.txt"

predictkcat_protein_sequences = project_root / "data_preprocessing" / "utils" / "predictkcat" / "fetch_ec_protein_sequences_only.py"
predictkcat_prepare = project_root / "data_preprocessing" / "utils" / "predictkcat" / "predictkcat_prepare.py"
predictkcat_predict = project_root / "data_preprocessing" / "utils" / "predictkcat" / "predictkcat_predict.py"
ecs_with_kcat_output_path = project_root / "data_preprocessing" / "utils" / "predictkcat" / "ecs.txt"

predictkcat_DLKcat = project_root / "data_preprocessing" / "utils" / "predictkcat" / "DLKcat"

mapping_ec_number_old_new = project_root / "data" / "raw" / "mapping_ec_numbers" / "mapping_ec_number_old_new.csv"
merged_overview = project_root / "data" / "processed" / "overview" / "overview.csv"

# ====== HYBRID MODEL OUTPUTS ======
reaction_index_csv               = project_root / "data" / "processed" / "overview" / "reaction_index.csv"
kcat_aggregate_csv               = project_root / "data" / "processed" / "overview" / "kcat_aggregate.csv"
hybrid_weights_csv               = project_root / "data" / "processed" / "overview" / "reaction_weights.csv"
ecfba_fluxes_csv                 = project_root / "data" / "processed" / "overview" / "ecfba_fluxes.csv"
active_pathways_json             = project_root / "data" / "processed" / "overview" / "predicted_pathways.json"

# ====== Model options ======
# Total enzyme budget (arbitrary units). If you have proteomics, scale appropriately.
E_TOTAL = 1.0

# Normalize & clamp for stability when converting kcat -> weights
KCAT_MIN = 1e-6     # s^-1 floor for numerical stability
KCAT_MAX = 1e6      # cap to suppress outliers

# Currency metabolite aliases
CURRENCY_ALIASES = {
    "H+": {"H+", "[H+]", "PROTON"},
    "H2O": {"H2O", "water"},
    "NAD+": {"NAD+", "nad+"},
    "NADH": {"NADH", "nadh"},
    "NADP+": {"NADP+", "nadp+"},
    "NADPH": {"NADPH", "nadph"},
    "ATP": {"ATP", "atp"},
    "ADP": {"ADP", "adp"},
    "Pi": {"phosphate", "Pi", "inorganic phosphate", "orthophosphate"},
}

# Default seed/target lists (you can edit when running)
DEFAULT_SEEDS = []   # e.g., ["glucose", "NAD+", "H2O"]
DEFAULT_TARGETS = [] # e.g., ["ethanol"]