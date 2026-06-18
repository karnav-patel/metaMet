# metaMet/data_preprocessing/hybrid/aggregate_kcat.py
from __future__ import annotations
import math
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Tuple


import data_preprocessing.config.config as config
from .utils_io import safe_read_csv

def _to_float(x):
    try:
        if x in (None, "", "NA", "nan"):
            return np.nan
        return float(x)
    except Exception:
        return np.nan

def _cap_kcat(x: float) -> float:
    if isinstance(x, float):
        if x <= 0 or math.isnan(x):
            return np.nan
        return float(min(max(x, config.KCAT_MIN), config.KCAT_MAX))
    return np.nan

def load_measured_kcats() -> pd.DataFrame:
    """
    Load measured kcat rows from BRENDA, MetaCyc, SABIO into a unified table:
    columns: ['ec_number', 'substrate', 'kcat', 'source']
    """
    rows = []

    # ---- BRENDA ----
    if Path(config.brenda_kcat_csv).exists():
        b = safe_read_csv(config.brenda_kcat_csv)
        for _, r in b.iterrows():
            ec = str(r.get("ec_number", "")).strip()
            sub = str(r.get("substrate", "")).strip()
            kc = _cap_kcat(_to_float(r.get("kcat")))
            if ec and not math.isnan(kc):
                rows.append((ec, sub, kc, "BRENDA"))

    # ---- MetaCyc ----
    if Path(config.metacyc_kcat_csv).exists():
        m = safe_read_csv(config.metacyc_kcat_csv)
        for _, r in m.iterrows():
            ec = str(r.get("ec_number", "")).strip()
            sub = str(r.get("substrate", "")).strip()
            kc = _cap_kcat(_to_float(r.get("kcat")))
            if ec and not math.isnan(kc):
                rows.append((ec, sub, kc, "MetaCyc"))

    # ---- SABIO ----
    if Path(config.sabio_kcat_csv).exists():
        s = safe_read_csv(config.sabio_kcat_csv, sep="\t")
        # SABIO columns seen: ECNumber, Substrate, Value, Unit
        for _, r in s.iterrows():
            ec = str(r.get("ECNumber", "")).strip()
            sub = str(r.get("Substrate", "")).strip()
            kc = _cap_kcat(_to_float(r.get("Value")))
            unit = str(r.get("Unit", "")).lower()
            # Convert if needed (we assume Value is in s^-1 as per your example)
            if ec and not math.isnan(kc):
                rows.append((ec, sub, kc, "SABIO"))

    df = pd.DataFrame(rows, columns=["ec_number", "substrate", "kcat", "source"])
    return df

def apply_ec_renumbering(df: pd.DataFrame, ec_col: str = "ec_number") -> pd.DataFrame:
    """
    Map old ECs to new ECs. If an old EC maps to multiple new ECs, duplicate rows.
    """
    map_df = safe_read_csv(config.mapping_ec_number_old_new)
    map_df.columns = ["old_ec", "new_ec"]
    map_df = map_df.dropna()
    m = df.copy()
    m["ec_number"] = m[ec_col].astype(str)
    m["_keep"] = True

    # Split into rows that need mapping and those that don't
    need = m[m["ec_number"].isin(set(map_df["old_ec"]))]
    keep = m[~m["ec_number"].isin(set(map_df["old_ec"]))]

    # Expand mappings
    expanded = need.merge(map_df, left_on="ec_number", right_on="old_ec", how="left")
    expanded["ec_number"] = expanded["new_ec"]
    expanded = expanded.drop(columns=["old_ec", "new_ec", "_keep"])

    out = pd.concat([keep.drop(columns=["_keep"]), expanded], ignore_index=True)
    out = out.dropna(subset=["ec_number"])
    return out

def aggregate_measured(df_meas: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse by ec_number -> median kcat and counts.
    """
    if df_meas.empty:
        return pd.DataFrame(columns=["ec_number", "kcat_measured_median", "n_measured"])
    g = df_meas.groupby("ec_number")["kcat"].agg(["median", "count"]).reset_index()
    g.columns = ["ec_number", "kcat_measured_median", "n_measured"]
    return g

def load_predicted_kcats() -> pd.DataFrame:
    """
    Join predictkcat input_with_ec (row order) with output.tsv (same row order),
    then aggregate by EC -> median predicted kcat and counts.
    """
    if not Path(config.predictkcat_input_with_ec).exists() or not Path(config.predictkcat_output_tsv).exists():
        return pd.DataFrame(columns=["ec_number", "kcat_pred_median", "n_pred"])

    inp = safe_read_csv(config.predictkcat_input_with_ec)
    outp = safe_read_csv(config.predictkcat_output_tsv, sep="\t")

    # Expect same ordering; if not, align by index length
    if len(inp) != len(outp):
        # Do a cautious inner join by (Substrate Name, Substrate SMILES, Protein Sequence)
        join_cols = ["Substrate Name", "Substrate SMILES", "Protein Sequence"]
        common = [c for c in join_cols if c in inp.columns and c in outp.columns]
        if common:
            merged = inp.merge(outp, on=common, how="inner")
        else:
            # last resort: truncate to min length
            n = min(len(inp), len(outp))
            inp = inp.iloc[:n].copy()
            outp = outp.iloc[:n].copy()
            merged = pd.concat([inp.reset_index(drop=True), outp.reset_index(drop=True)], axis=1)
    else:
        merged = pd.concat([inp.reset_index(drop=True), outp.reset_index(drop=True)], axis=1)

    # Output has a column like 'Kcat value (1/s)'
    kc_col = None
    for c in merged.columns:
        if "kcat" in c.lower() and "/s" in c.lower():
            kc_col = c
            break
    if kc_col is None:
        # try generic
        kc_col = "Kcat value (1/s)" if "Kcat value (1/s)" in merged.columns else merged.columns[-1]

    merged["kcat_pred"] = merged[kc_col].apply(_to_float).apply(_cap_kcat)
    merged = merged.dropna(subset=["ec_number", "kcat_pred"])
    g = merged.groupby("ec_number")["kcat_pred"].agg(["median", "count"]).reset_index()
    g.columns = ["ec_number", "kcat_pred_median", "n_pred"]
    return g

def build_kcat_aggregate() -> pd.DataFrame:
    """
    Combine measured and predicted kcat per EC; choose a final per-EC kcat:
      - prefer measured median if available; otherwise predicted median.
    Write:
      - kcat_aggregate_csv
      - merged_ec_numbers_with_kcat (txt)
      - predictkcat_missing_kcat_ids_output_path (txt)
      - merged_overview (csv, includes counts)
    """
    # Load + map measured
    meas = load_measured_kcats()
    if not meas.empty:
        meas = apply_ec_renumbering(meas)
    meas_agg = aggregate_measured(meas)

    # Load + map predicted
    pred = load_predicted_kcats()
    if not pred.empty:
        pred = apply_ec_renumbering(pred, ec_col="ec_number")
    pred_agg = pred

    # Merge
    all_ec = pd.merge(meas_agg, pred_agg, on="ec_number", how="outer")
    # Choose final
    all_ec["kcat_final"] = all_ec["kcat_measured_median"].where(
        ~all_ec["kcat_measured_median"].isna(),
        all_ec["kcat_pred_median"]
    )
    all_ec["has_measured"] = ~all_ec["kcat_measured_median"].isna()
    all_ec["has_pred"] = ~all_ec["kcat_pred_median"].isna()

    # Save kcat aggregate
    Path(config.kcat_aggregate_csv).parent.mkdir(parents=True, exist_ok=True)
    all_ec.to_csv(config.kcat_aggregate_csv, index=False)

    # Load merged EC ids universe (all)
    all_ids = []
    for p in [
        config.brenda_extracted_ids_output_path,
        config.kegg_extracted_ids_output_path,
        config.metacyc_extracted_ids_output_path,
    ]:
        if Path(p).exists():
            ids = [s.strip() for s in Path(p).read_text().splitlines() if s.strip()]
            all_ids.extend(ids)
    all_ids = sorted(set(all_ids))
    Path(config.merged_extracted_ids_output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(config.merged_extracted_ids_output_path).write_text("\n".join(all_ids) + "\n")

    # ECs that have any kcat
    ec_with_kcat = sorted(set(all_ec[~all_ec["kcat_final"].isna()]["ec_number"]))
    Path(config.merged_ec_numbers_with_kcat).write_text("\n".join(ec_with_kcat) + "\n")
    Path(config.ecs_with_kcat_output_path).write_text("\n".join(ec_with_kcat) + "\n")

    # ECs missing kcat
    missing = [ec for ec in all_ids if ec not in set(ec_with_kcat)]
    Path(config.predictkcat_missing_kcat_ids_output_path).write_text("\n".join(missing) + "\n")

    # High-level merged overview
    overview = all_ec.rename(columns={
        "kcat_measured_median": "kcat_measured_median_s^-1",
        "kcat_pred_median": "kcat_pred_median_s^-1",
        "kcat_final": "kcat_selected_s^-1"
    })
    overview.to_csv(config.merged_overview, index=False)
    return overview
