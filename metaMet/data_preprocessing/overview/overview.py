#!/usr/bin/env python3
"""
Build a high-level EC summary CSV combining:
  • Measured kcats from BRENDA, MetaCyc, Sabio
  • Reaction entries from BRENDA, MetaCyc, KEGG
  • Model-predicted kcats (aligned by row number to input_with_ec)

Rules:
  – Count a reaction only if it has both educts & products.
  – Count a kcat only if at least one valid numeric measurement exists.
  – Handle EC renumbering via mapping (old → new).

Output: config.merged_overview
"""
from __future__ import annotations

import sys
import re
from pathlib import Path
import pandas as pd


# --- Locate project root and import config -----------------------------------
def get_project_root() -> Path:
    current_file_path = Path(__file__).resolve()
    for parent in current_file_path.parents:
        if parent.name == "metaMet":
            return parent
    raise RuntimeError("Project root folder 'metaMet' not found.")

project_root = get_project_root()
data_preproc_dir = project_root / "data_preprocessing"
sys.path.insert(0, str(data_preproc_dir))

try:
    import config.config as config
except Exception as e:
    raise RuntimeError(
        "Could not import config.config. Ensure it is available under "
        f"{data_preproc_dir}.\nOriginal error: {e}"
    )


# === CONFIG (all from config.config; no fallbacks) ============================
try:
    MASTER_EC_FILE = Path(config.merged_extracted_ids_output_path)
    MAPPING_FILE   = Path(config.mapping_ec_number_old_new)

    BRENDA_KCAT_FILE = Path(config.brenda_kcat_csv)
    BRENDA_RXN_FILE  = Path(config.brenda_reaction_csv)

    METACYC_KCAT_FILE = Path(config.metacyc_kcat_csv)
    METACYC_RXN_FILE  = Path(config.metacyc_reaction_csv)

    KEGG_RXN_FILE = Path(config.kegg_reaction_csv)

    SABIO_KCAT_FILE = Path(config.sabio_kcat_csv)

    PREDICTKCAT_INPUT_TSV     = Path(config.predictkcat_input_tsv)           # not directly used; kept for sanity
    PREDICTKCAT_INPUT_WITH_EC = Path(config.predictkcat_input_with_ec)
    PREDICTKCAT_OUTPUT_TSV    = Path(config.predictkcat_output_tsv)

    OUTPUT_FILE = Path(config.merged_overview)
except AttributeError as e:
    missing = str(e).split("'")[-2]
    raise RuntimeError(f"Missing required config field: config.{missing}")


# === Helpers =================================================================
def require_exists(p: Path, label: str):
    if not p.exists():
        raise RuntimeError(f"{label} not found at: {p}")

def load_mapping(path: Path):
    """
    Mapping CSV with two columns: 'old ec number' and 'new number' (any case/whitespace),
    or exactly two columns -> [old_ec, new_ec].
    Returns: (old2new: dict[str,set[str]], new2old: dict[str,set[str]])
    """
    df = pd.read_csv(path, dtype=str).fillna("")
    cols = [c.strip().lower() for c in df.columns]
    if "old ec number" in cols and "new number" in cols:
        df = df.rename(columns={
            df.columns[cols.index("old ec number")]: "old_ec",
            df.columns[cols.index("new number")]:    "new_ec"
        })
    elif len(df.columns) == 2:
        df.columns = ["old_ec", "new_ec"]
    else:
        raise RuntimeError(
            f"Unexpected columns in mapping file {path}. "
            "Provide either ['old ec number','new number'] or exactly 2 columns."
        )
    df["old_ec"] = df["old_ec"].str.strip()
    df["new_ec"] = df["new_ec"].str.strip()

    old2new = df.groupby("old_ec")["new_ec"].apply(lambda s: set(s) - {""}).to_dict()
    new2old: dict[str, set[str]] = {}
    for o, news in old2new.items():
        for n in news:
            new2old.setdefault(n, set()).add(o)
    return old2new, new2old

def get_synonyms(ec: str, old2new: dict[str,set[str]], new2old: dict[str,set[str]]) -> set[str]:
    s = {ec}
    s |= old2new.get(ec, set())
    s |= new2old.get(ec, set())
    return s

def load_master_ecs(path: Path) -> set[str]:
    require_exists(path, "Master EC list")
    # handle either txt (one EC per line) or csv with EC in first column
    try:
        df = pd.read_csv(path, header=None, dtype=str)
    except Exception:
        df = pd.read_csv(path, header=0, dtype=str)
    ec_series = df.iloc[:, 0].astype(str).str.strip()
    return set(e for e in ec_series if e and e.lower() != "ec")


def load_brenda_kcat() -> set[str]:
    require_exists(BRENDA_KCAT_FILE, "BRENDA kcat CSV")
    df = pd.read_csv(BRENDA_KCAT_FILE, dtype=str).fillna("")
    if "ec_number" not in df.columns:
        return set()
    if "kcat" in df.columns:
        knum = pd.to_numeric(df["kcat"], errors="coerce")
        return set(df.loc[~knum.isna(), "ec_number"].astype(str).str.strip())
    kcat_cols = [c for c in df.columns if "kcat" in c.lower()]
    if not kcat_cols:
        return set()
    mask = False
    for c in kcat_cols:
        mask = mask | pd.to_numeric(df[c], errors="coerce").notna()
    return set(df.loc[mask, "ec_number"].astype(str).str.strip())

def load_metacyc_kcat() -> set[str]:
    require_exists(METACYC_KCAT_FILE, "MetaCyc kcat CSV")
    df = pd.read_csv(METACYC_KCAT_FILE, dtype=str, engine="python", quotechar='"').fillna("")
    if "ec_number" not in df.columns:
        return set()
    kcat_cols = [c for c in df.columns if "kcat" in c.lower()]
    if not kcat_cols:
        return set()
    mask = False
    for c in kcat_cols:
        mask = mask | pd.to_numeric(df[c], errors="coerce").notna()
    ec_pattern = re.compile(r"^\d+(\.\d+){3}$")
    df["ec_number"] = df["ec_number"].astype(str).str.strip()
    valid_ec = df["ec_number"].apply(lambda x: bool(ec_pattern.match(x)))
    good = df.loc[valid_ec & mask, "ec_number"]
    return set(good)

def load_sabio_kcat() -> set[str]:
    require_exists(SABIO_KCAT_FILE, "SABIO kcat TSV")
    df = pd.read_csv(SABIO_KCAT_FILE, sep="\t", dtype=str).fillna("")
    if "Type" not in df.columns or "ECNumber" not in df.columns:
        return set()
    df = df[df["Type"].str.lower() == "kcat"]
    return set(df["ECNumber"].astype(str).str.strip())

def _load_rxn_file(path: Path, label: str, ec_col="ec_number", ed_col="educts", pr_col="products") -> set[str]:
    require_exists(path, label)
    df = pd.read_csv(path, dtype=str).fillna("")
    if not {ec_col, ed_col, pr_col}.issubset(df.columns):
        return set()
    ok = df[ed_col].astype(bool) & df[pr_col].astype(bool)
    return set(df.loc[ok, ec_col].astype(str).str.strip())

def load_brenda_reaction() -> set[str]:
    return _load_rxn_file(BRENDA_RXN_FILE, "BRENDA reactions CSV")

def load_metacyc_reaction() -> set[str]:
    return _load_rxn_file(METACYC_RXN_FILE, "MetaCyc reactions CSV")

def load_kegg_reaction() -> set[str]:
    return _load_rxn_file(KEGG_RXN_FILE, "KEGG reactions CSV")


def load_predicted_ecs_by_row() -> set[str]:
    """
    Determine which ECs have a model-predicted kcat by aligning
    predictkcat_output_tsv rows to predictkcat_input_with_ec rows by index.
    """
    require_exists(PREDICTKCAT_INPUT_WITH_EC, "predictkcat input_with_ec CSV")
    require_exists(PREDICTKCAT_OUTPUT_TSV, "predictkcat output TSV")

    inp = pd.read_csv(PREDICTKCAT_INPUT_WITH_EC, dtype=str).fillna("")
    if "ec_number" not in inp.columns:
        raise RuntimeError(
            f"'ec_number' column not found in {PREDICTKCAT_INPUT_WITH_EC}."
        )

    out = pd.read_csv(PREDICTKCAT_OUTPUT_TSV, sep="\t", dtype=str).fillna("")

    # Prefer columns containing 'kcat', else allow typical prediction-like names
    kcat_cols = [c for c in out.columns if "kcat" in c.lower()]
    if not kcat_cols:
        kcat_cols = [c for c in out.columns if any(t in c.lower() for t in ["pred", "rate", "value"])]
    if not kcat_cols:
        raise RuntimeError(
            f"No predicted-kcat-like column found in {PREDICTKCAT_OUTPUT_TSV} "
            "(expected a column containing 'kcat', 'pred', 'rate', or 'value')."
        )

    # Any numeric in these columns means we have a prediction for that row
    pred_mask = False
    for c in kcat_cols:
        pred_mask = pred_mask | pd.to_numeric(out[c], errors="coerce").notna()

    if len(inp) != len(out):
        raise RuntimeError(
            f"Row mismatch: input_with_ec has {len(inp)} rows, output.tsv has {len(out)} rows. "
            "They must align by row index."
        )

    ecs_with_pred = set(inp.loc[pred_mask, "ec_number"].astype(str).str.strip())
    ecs_with_pred.discard("")
    return ecs_with_pred


def build_summary(master: set[str],
                  old2new: dict[str,set[str]],
                  new2old: dict[str,set[str]],
                  bk: set[str], mk: set[str], sk: set[str],
                  br: set[str], mr: set[str], kr: set[str],
                  pe: set[str]) -> pd.DataFrame:
    all_ecs = set(master) | bk | mk | sk | br | mr | kr | pe
    for e in list(all_ecs):
        all_ecs |= old2new.get(e, set())

    rows = []
    for ec in sorted(all_ecs):
        syn = get_synonyms(ec, old2new, new2old)

        ksrc = []
        if syn & bk: ksrc.append("BRENDA")
        if syn & mk: ksrc.append("MetaCyc")
        if syn & sk: ksrc.append("Sabio")

        rsrc = []
        if syn & br: rsrc.append("BRENDA")
        if syn & mr: rsrc.append("MetaCyc")
        if syn & kr: rsrc.append("KEGG")

        rows.append({
            "all_EC":          ec,
            "newEC":           "|".join(sorted(old2new.get(ec, {ec}))),
            "existing_kcat":   "Yes" if ksrc else "No",
            "kcat_source":     "|".join(ksrc),
            "having_reaction": "Yes" if rsrc else "No",
            "reaction_source": "|".join(rsrc),
            "predicted_kcat":  "Yes" if (syn & pe) else "No",
        })
    return pd.DataFrame(rows)


def main():
    # Load inputs
    master_ecs = load_master_ecs(MASTER_EC_FILE)
    old2new, new2old = load_mapping(MAPPING_FILE)

    br_k = load_brenda_kcat()
    mc_k = load_metacyc_kcat()
    sb_k = load_sabio_kcat()

    br_r = load_brenda_reaction()
    mc_r = load_metacyc_reaction()
    kg_r = load_kegg_reaction()

    pred_ec = load_predicted_ecs_by_row()

    # Build and write summary
    summary = build_summary(
        master_ecs, old2new, new2old,
        br_k, mc_k, sb_k,
        br_r, mc_r, kg_r,
        pred_ec
    )

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUTPUT_FILE, index=False)
    print(f"Wrote summary to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
