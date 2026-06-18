# metaMet/data_preprocessing/hybrid/utils_io.py
from __future__ import annotations
import ast
import math
import pandas as pd
from typing import Any, List

def parse_list_cell(cell: Any) -> List[str]:
    """
    Parse a cell that looks like "['A','B']" into a Python list.
    Falls back to [str(cell)] if parsing fails or cell is empty/NaN.
    """
    if cell is None or (isinstance(cell, float) and math.isnan(cell)):
        return []
    s = str(cell).strip()
    if s == "":
        return []
    try:
        val = ast.literal_eval(s)
        if isinstance(val, list):
            return [str(x).strip() for x in val]
        return [str(val)]
    except Exception:
        return [s]

def safe_read_csv(path, **kwargs) -> pd.DataFrame:
    df = pd.read_csv(path, **kwargs)
    # Normalize column names
    df.columns = [c.strip() for c in df.columns]
    return df

def normalize_name(name: str) -> str:
    """
    Basic name normalization: lower, strip, collapse spaces.
    Leave SMILES alone (handled separately).
    """
    s = (name or "").strip().lower()
    s = " ".join(s.split())
    return s
