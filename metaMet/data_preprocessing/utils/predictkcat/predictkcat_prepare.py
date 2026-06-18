# predictkcat_prepare_from_seq_csv.py
import os
import re
import sys
import csv
import ast
import time
import requests
import pubchempy as pcp
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

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

def get_ecs_missing_kcat():
    brenda_file = config.brenda_kcat_csv
    metacyc_file = config.metacyc_kcat_csv
    sabio_file  = config.sabio_kcat_csv
    merged_file = config.merged_extracted_ids_output_path

    all_ecs = set()
    with merged_file.open("r", encoding="utf-8") as f:
        for line in f:
            ec = line.strip()
            if ec:
                all_ecs.add(ec)

    ecs_with_kcat = set()

    def iter_rows(path: Path, delimiter=","):
        with path.open("r", encoding="utf-8", newline="") as f:
            yield from csv.DictReader(f, delimiter=delimiter)

    for row in iter_rows(brenda_file, delimiter=","):
        ec   = row.get("ec_number", "").strip()
        kcat = row.get("kcat", "").strip()
        if ec and kcat:
            ecs_with_kcat.add(ec)

    for row in iter_rows(metacyc_file, delimiter=","):
        ec   = row.get("ec_number", "").strip()
        kcat = row.get("kcat", "").strip()
        if ec and kcat:
            ecs_with_kcat.add(ec)

    for row in iter_rows(sabio_file, delimiter="\t"):
        if row.get("Type", "").strip().lower() == "kcat":
            ec = row.get("ECNumber", "").strip()
            if ec:
                ecs_with_kcat.add(ec)

    missing_ecs = sorted(all_ecs - ecs_with_kcat)
    return missing_ecs, sorted(ecs_with_kcat)

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

def _get_txt(url):
    r = requests.get(url, timeout=5)
    r.raise_for_status()
    return r.text.strip()

def fetch_smiles(name):
    name_q = requests.utils.quote(name)
    try:
        return _get_txt(f"https://opsin.ch.cam.ac.uk/opsin/{name_q}.smiles")
    except: pass
    try:
        return _get_txt(f"{PUBCHEM}/compound/name/{name_q}/property/CanonicalSMILES/TXT")
    except: pass
    try:
        return _get_txt(f"http://cactus.nci.nih.gov/chemical/structure/{name_q}/smiles")
    except: pass
    try:
        comps = pcp.get_compounds(name, "name")
        if comps: return comps[0].canonical_smiles or ""
    except: pass
    if name.upper() in ("H+","[H+]"): return "[H+]"
    return ""

def _clean_fields(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    # Substrate names: remove tabs/newlines/quotes, trim
    df["substrate"] = (
        df["substrate"].astype(str)
        .str.replace(r'[\t\r\n"]+', " ", regex=True)
        .str.strip()
    )
    # SMILES: remove all whitespace and quotes (SMILES must be 1 token)
    df["smiles"] = (
        df["smiles"].astype(str)
        .str.replace(r'["\s]+', "", regex=True)
        .str.strip()
    )
    # Protein sequences: keep only letters, uppercase, no spaces
    df["protein_sequence"] = (
        df["protein_sequence"].astype(str)
        .str.replace(r'[^A-Za-z]', "", regex=True)
        .str.upper()
        .str.strip()
    )
    return df

def build_input_and_output():
    # 1) Take last 5 missing ECs (but write full lists to disk)
    missing_list, ecs_with_kcat = get_ecs_missing_kcat()
    if not missing_list:
        print("→ No missing ECs to process. Exiting.")
        # Still write the with-kcat file for auditing
        with open(config.merged_ec_numbers_with_kcat, "w", encoding="utf-8") as f:
            for ec in ecs_with_kcat:
                f.write(f"{ec}\n")
        sys.exit(0)

    # Write the full lists
    missing_txt = config.predictkcat_missing_kcat_ids_output_path
    with missing_txt.open("w", encoding="utf-8") as f:
        for ec in missing_list:
            f.write(f"{ec}\n")

    with open(config.merged_ec_numbers_with_kcat, "w", encoding="utf-8") as f:
        for ec in ecs_with_kcat:
            f.write(f"{ec}\n")

    # 2) Load Metacyc smiles-enriched file
    smiles_path = config.metacyc_reaction_csv_with_smiles
    df_smiles = pd.read_csv(smiles_path, dtype=str).fillna("")
    df_smiles["ec_clean"] = df_smiles["ec_number"].str.strip().str.lower()

    # 3) Load BRENDA, KEGG, MetaCyc reaction tables
    brenda_df = pd.read_csv(config.brenda_reaction_csv, dtype=str).fillna("")
    brenda_df["ec_clean"] = brenda_df["ec_number"].str.strip().str.lower()

    kegg_df = pd.read_csv(config.kegg_reaction_csv, dtype=str).fillna("")
    kegg_df["ec_clean"] = kegg_df["ec_number"].str.strip().str.lower()

    metacyc_df = pd.read_csv(config.metacyc_reaction_csv, dtype=str).fillna("")
    metacyc_df["ec_clean"] = metacyc_df["ec_number"].str.strip().str.lower()

    rows = []
    seen = set()

    for ec in missing_list:
        ec_clean = ec.strip().lower()

        # 3a) Prefer MetaCyc with SMILES (pad-zip so missing SMILES still produce rows)
        meta_hits = df_smiles[df_smiles["ec_clean"] == ec_clean]
        added_any = False
        if not meta_hits.empty:
            for _, r in meta_hits.iterrows():
                try:
                    eds = ast.literal_eval(r["educts"])
                    es  = ast.literal_eval(r["educts_smiles"])
                    pds = ast.literal_eval(r["products"])
                    ps  = ast.literal_eval(r["products_smiles"])
                except Exception:
                    continue

                for name_list, smile_list in ((eds, es), (pds, ps)):
                    max_len = max(len(name_list), len(smile_list))
                    for i in range(max_len):
                        name = name_list[i].strip() if i < len(name_list) else ""
                        smi  = smile_list[i].strip() if i < len(smile_list) else ""
                        if not name:  # skip truly empty names
                            continue
                        key = (ec_clean, name)
                        if key not in seen:
                            seen.add(key)
                            rows.append({
                                "ec_number": ec_clean,
                                "substrate": name,
                                "smiles":    smi  # might be "", API fill happens later
                            })
                            added_any = True

            # Only skip fallbacks if we actually captured something from MetaCyc
            if added_any:
                continue


        # 3b) Fallback to BRENDA → KEGG → MetaCyc (no SMILES)
        for df in (brenda_df, kegg_df, metacyc_df):
            hits = df[df["ec_clean"] == ec_clean]
            for _, r in hits.iterrows():
                try:
                    eds = ast.literal_eval(r["educts"])
                    pds = ast.literal_eval(r["products"])
                except:
                    continue
                for name in eds + pds:
                    key = (ec_clean, name)
                    if key not in seen:
                        seen.add(key)
                        rows.append({
                            "ec_number": ec_clean,
                            "substrate": name.strip(),
                            "smiles":    ""  # fetch later
                        })

    # 4) Fetch missing SMILES
    for row in rows:
        if not row["smiles"]:
            row["smiles"] = fetch_smiles(row["substrate"])
            time.sleep(0.2)

    # Clean SMILES and validate
    for row in rows:
        smi = row["smiles"].strip()
        smi = re.sub(r"\[R\d*\]", "*", smi)
        mol = Chem.MolFromSmiles(smi)
        row["smiles"] = smi if mol else ""

    # 5) Read protein sequences from CSV (NO online fetching here)
    seq_df = pd.read_csv(config.predictkcat_with_ec_protein_sequences, dtype=str).fillna("")
    seq_df["ec_number"] = seq_df["ec_number"].str.strip().str.lower()
    seq_df["protein_sequence"] = seq_df["protein_sequence"].str.strip()
    seq_map = dict(zip(seq_df["ec_number"], seq_df["protein_sequence"]))

    # 6) Assemble DataFrame
    df_out = pd.DataFrame(rows)
    df_out["protein_sequence"] = df_out["ec_number"].map(seq_map)

    # 7) Filter complete rows
    complete_df = df_out[
        (df_out["substrate"] != "") &
        (df_out["smiles"] != "") &
        (df_out["protein_sequence"].notna()) &
        (df_out["protein_sequence"] != "")
    ]

    complete_df = _clean_fields(complete_df)

    # 8) Write outputs
    input_tsv  = config.predictkcat_input_tsv
    output_csv = config.predictkcat_input_with_ec
    failed_txt = config.predictkcat_failed_kcat_ids_output_path

    with open(input_tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n",
                    quoting=csv.QUOTE_NONE, escapechar="\\")
        w.writerow(["Substrate Name", "Substrate SMILES", "Protein Sequence"])
        for _, r in complete_df.iterrows():
            w.writerow([r["substrate"], r["smiles"], r["protein_sequence"]])

    complete_df.rename(columns={
        "ec_number":        "ec_number",
        "substrate":        "Substrate Name",
        "smiles":           "Substrate SMILES",
        "protein_sequence": "Protein Sequence"
    })[["ec_number", "Substrate Name", "Substrate SMILES", "Protein Sequence"]].to_csv(
        output_csv, index=False
    )

    # 9) Record ECs that failed to produce any complete row
    failed_ecs = sorted(set(df_out["ec_number"]) - set(complete_df["ec_number"]))
    with failed_txt.open("w") as f:
        for ec in failed_ecs:
            f.write(f"{ec}\n")

if __name__ == "__main__":
    build_input_and_output()
