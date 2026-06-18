# metaMet/data_preprocessing/hybrid/build_reaction_index.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd


import data_preprocessing.config.config as config
from .utils_io import safe_read_csv, parse_list_cell, normalize_name

def _canonicalize_metabolite(name: str, smiles: str | None, alias_map: Dict[str, str]) -> str:
    """
    Prefer SMILES identifier if provided; else normalized name.
    Apply alias unification for currency metabolites.
    """
    if smiles and isinstance(smiles, str) and smiles.strip().startswith(("C", "O", "N", "S", "[")):
        return f"SMI::{smiles.strip()}"
    n = normalize_name(name)
    # unify currencies to canonical keys
    for canon, alset in config.CURRENCY_ALIASES.items():
        if n in {a.lower() for a in alset} or n == canon.lower():
            return f"CUR::{canon}"
    return f"TXT::{n}"

def load_metacyc_smiles_map() -> Dict[str, str]:
    """
    Build a quick map of 'raw metabolite name' -> SMILES from the MetaCyc file if present.
    """
    m = {}
    if Path(config.metacyc_reaction_csv_with_smiles).exists():
        df = safe_read_csv(config.metacyc_reaction_csv_with_smiles)
        # columns: educts, products, educts_smiles, products_smiles
        for _, r in df.iterrows():
            ed = parse_list_cell(r.get("educts"))
            pr = parse_list_cell(r.get("products"))
            es = parse_list_cell(r.get("educts_smiles"))
            ps = parse_list_cell(r.get("products_smiles"))
            for n, s in list(zip(ed, es)) + list(zip(pr, ps)):
                if n and s:
                    m[str(n).strip()] = str(s).strip()
    return m

def _pack_rows(rows: List[Tuple]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=[
        "reaction_id","source","ec_number","educts","products"
    ])

def build_reaction_index() -> pd.DataFrame:
    """
    Read reactions from BRENDA/KEGG/MetaCyc into a unified CSV.
    Using unit stoichiometry (1 per metabolite) as a scaffold.
    """
    meta_smiles = load_metacyc_smiles_map()
    rows = []

    # --- BRENDA ---
    if Path(config.brenda_reaction_csv).exists():
        df = safe_read_csv(config.brenda_reaction_csv)
        for i, r in df.iterrows():
            ec = str(r.get("ec_number","")).strip()
            ed = parse_list_cell(r.get("educts"))
            pr = parse_list_cell(r.get("products"))
            if not ec or not ed or not pr:
                continue
            rid = f"BRENDA::{ec}::{i}"
            rows.append((rid, "BRENDA", ec, ed, pr))

    # --- KEGG ---
    if Path(config.kegg_reaction_csv).exists():
        df = safe_read_csv(config.kegg_reaction_csv)
        for i, r in df.iterrows():
            ec = str(r.get("ec_number","")).strip()
            ed = parse_list_cell(r.get("educts"))
            pr = parse_list_cell(r.get("products"))
            if not ec or not ed or not pr:
                continue
            rid = f"KEGG::{ec}::{i}"
            rows.append((rid, "KEGG", ec, ed, pr))

    # --- MetaCyc ---
    if Path(config.metacyc_reaction_csv).exists():
        df = safe_read_csv(config.metacyc_reaction_csv)
        for i, r in df.iterrows():
            ec = str(r.get("ec_number","")).strip()
            ed = parse_list_cell(r.get("educts"))
            pr = parse_list_cell(r.get("products"))
            if not ec or not ed or not pr:
                continue
            rid = f"MetaCyc::{ec}::{r.get('ReactionID', i)}"
            rows.append((rid, "MetaCyc", ec, ed, pr))

    rxn_df = _pack_rows(rows)

    # Canonicalize metabolite identifiers (prefer SMILES for MetaCyc if available)
    canon_ed, canon_pr = [], []
    for _, row in rxn_df.iterrows():
        ed = []
        for n in row["educts"]:
            smi = meta_smiles.get(n)
            ed.append(_canonicalize_metabolite(n, smi, config.CURRENCY_ALIASES))
        pr = []
        for n in row["products"]:
            smi = meta_smiles.get(n)
            pr.append(_canonicalize_metabolite(n, smi, config.CURRENCY_ALIASES))
        canon_ed.append(ed); canon_pr.append(pr)

    rxn_df = rxn_df.assign(educts=canon_ed, products=canon_pr)

    Path(config.reaction_index_csv).parent.mkdir(parents=True, exist_ok=True)
    rxn_df.to_csv(config.reaction_index_csv, index=False)
    return rxn_df
