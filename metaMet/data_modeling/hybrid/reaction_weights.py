# data_modeling/hybrid/reaction_weights.py
from __future__ import annotations
import math
import numpy as np
import pandas as pd
from pathlib import Path
import data_preprocessing.config.config as config
from .utils_io import safe_read_csv

def _find_col(df: pd.DataFrame, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def make_reaction_weights(kcat_aggregate_csv: str, reaction_index_csv: str) -> pd.DataFrame:
    rxn = safe_read_csv(reaction_index_csv)
    kc  = safe_read_csv(kcat_aggregate_csv)

    # tolerate both pre- and post-rename schemas
    kcat_col = _find_col(kc, ["kcat_selected_s^-1", "kcat_final", "kcat_selected"])
    if kcat_col is None:
        raise KeyError(
            "No selected kcat column found in kcat aggregate. "
            f"Available columns: {list(kc.columns)}. "
            "Expected one of: 'kcat_selected_s^-1' (preferred) or 'kcat_final'/'kcat_selected'."
        )

    has_measured_col = "has_measured" if "has_measured" in kc.columns else None
    has_pred_col     = "has_pred"     if "has_pred"     in kc.columns else None

    cols = ["ec_number", kcat_col]
    if has_measured_col: cols.append(has_measured_col)
    if has_pred_col:     cols.append(has_pred_col)

    sub = kc[cols].copy().rename(columns={kcat_col: "kcat_selected_s^-1"})
    if "has_measured" not in sub.columns: sub["has_measured"] = False
    if "has_pred"     not in sub.columns: sub["has_pred"]     = False

    m = rxn.merge(sub, on="ec_number", how="left")

    def _logw(x):
        try:
            x = float(x)
            if not np.isfinite(x) or x <= 0: return np.nan
            x = min(max(x, config.KCAT_MIN), config.KCAT_MAX)
            return math.log10(x)
        except Exception:
            return np.nan

    m["log10_kcat"] = m["kcat_selected_s^-1"].apply(_logw)

    finite = m["log10_kcat"].replace([np.inf, -np.inf], np.nan).dropna()
    if len(finite) > 0:
        lo, hi = np.percentile(finite, [5, 95])
        span = max(hi - lo, 1e-6)
        def _norm(z):
            if z is None or (isinstance(z, float) and np.isnan(z)): return 0.05
            return float(min(max((z - lo) / span, 0.0), 1.0))
        m["weight"] = m["log10_kcat"].apply(_norm)
    else:
        m["weight"] = 0.05

    # small trust bonus for measured kcat
    m["weight"] = m["weight"] + 0.1 * m["has_measured"].fillna(False).astype(float)
    m["weight"] = m["weight"].clip(0, 1.5)

    Path(config.hybrid_weights_csv).parent.mkdir(parents=True, exist_ok=True)
    m.to_csv(config.hybrid_weights_csv, index=False)
    return m
