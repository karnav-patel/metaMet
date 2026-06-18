# metaMet/modeling/ec_fba.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from pathlib import Path

from scipy import sparse  # NEW

import data_preprocessing.config.config as config
from hybrid.utils_io import safe_read_csv

@dataclass
class EcFBAResult:
    reactions: List[str]
    metabolites: List[str]
    vf: np.ndarray
    vr: np.ndarray
    ef: np.ndarray
    er: np.ndarray

def build_stoichiometry(rxn_df: pd.DataFrame) -> Tuple[List[str], List[str], sparse.csr_matrix]:
    """
    Build sparse S (metabolites x reactions) using unit stoichiometry.
    """
    mets = set()
    for _, r in rxn_df.iterrows():
        for m in r["educts"]:
            mets.add(m)
        for m in r["products"]:
            mets.add(m)
    mets = sorted(mets)
    met_index = {m: i for i, m in enumerate(mets)}

    R = len(rxn_df)
    rows, cols, data = [], [], []
    for j, (_, r) in enumerate(rxn_df.iterrows()):
        for m in r["educts"]:
            rows.append(met_index[m]); cols.append(j); data.append(-1.0)
        for m in r["products"]:
            rows.append(met_index[m]); cols.append(j); data.append(+1.0)

    S = sparse.coo_matrix((data, (rows, cols)),
                          shape=(len(mets), R),
                          dtype=np.float64).tocsr()
    return mets, list(rxn_df["reaction_id"].tolist()), S

def _linprog_solve(c, Aeq, beq, Aub, bub, lb, ub):
    """
    Solve LP: minimize c^T x s.t. Aeq x = beq; Aub x <= bub; lb <= x <= ub
    using scipy.optimize.linprog (HiGHS). Supports sparse matrices.
    """
    from scipy.optimize import linprog
    bounds = list(zip(lb, ub))
    # Better to use None instead of np.inf for unbounded
    bounds = [(l, None if (u is None or np.isinf(u)) else u) for (l, u) in bounds]
    res = linprog(c, A_eq=Aeq, b_eq=beq, A_ub=Aub, b_ub=bub,
                  bounds=bounds, method="highs")
    return res

def solve_ecfba(
    rxn_df: pd.DataFrame,
    weights_df: pd.DataFrame,
    seeds: List[str] | None = None,
    targets: List[str] | None = None,
    E_total: float = None
) -> EcFBAResult:

    if E_total is None:
        E_total = float(config.E_TOTAL)

    m = rxn_df.merge(
        weights_df[["reaction_id","weight","kcat_selected_s^-1"]],
        on="reaction_id", how="left"
    )
    m["weight"] = m["weight"].fillna(0.05).astype(float)
    m["kcat_selected_s^-1"] = m["kcat_selected_s^-1"].fillna(config.KCAT_MIN).astype(float)

    mets, rxns, S = build_stoichiometry(m)

    R = len(rxns)
    VF = slice(0, R)
    VR = slice(R, 2*R)
    EF = slice(2*R, 3*R)
    ER = slice(3*R, 4*R)
    nvar = 4*R

    # Objective: maximize w·(vf+vr) == minimize -w·(vf+vr)
    w = m["weight"].to_numpy(dtype=np.float64)
    c = np.zeros(nvar, dtype=np.float64)
    c[VF] = -w
    c[VR] = -w

    # Aeq x = beq  -> [ S | -S | 0 | 0 ]
    Z_met_R = sparse.csr_matrix((len(mets), R), dtype=np.float64)
    Aeq = sparse.hstack([S, -S, Z_met_R, Z_met_R], format="csr")
    beq = np.zeros(len(mets), dtype=np.float64)

    # Aub x <= bub:
    # vf <= kcat*ef : [ I | 0 | -diag(kcat) | 0 ]
    # vr <= kcat*er : [ 0 | I | 0 | -diag(kcat) ]
    I = sparse.eye(R, dtype=np.float64, format="csr")
    O = sparse.csr_matrix((R, R), dtype=np.float64)
    kcat = m["kcat_selected_s^-1"].to_numpy(dtype=np.float64)
    K = sparse.diags(kcat, offsets=0, shape=(R, R), dtype=np.float64, format="csr")

    row_vf = sparse.hstack([I, O, -K, O], format="csr")
    row_vr = sparse.hstack([O, I, O, -K], format="csr")

    # Enzyme budget: sum(ef + er) <= E_total
    cols_budget = np.concatenate([np.arange(EF.start, EF.stop), np.arange(ER.start, ER.stop)])
    data_budget = np.ones(cols_budget.size, dtype=np.float64)
    rows_budget = np.zeros(cols_budget.size, dtype=np.int64)
    budget = sparse.coo_matrix((data_budget, (rows_budget, cols_budget)),
                               shape=(1, nvar), dtype=np.float64).tocsr()

    Aub = sparse.vstack([row_vf, row_vr, budget], format="csr")
    bub = np.concatenate([np.zeros(2*R, dtype=np.float64), np.array([float(E_total)], dtype=np.float64)])

    # Bounds: x >= 0 (no finite upper bounds)
    lb = np.zeros(nvar, dtype=np.float64)
    ub = np.full(nvar, None, dtype=object)  # None == +inf to linprog

    res = _linprog_solve(c, Aeq, beq, Aub, bub, lb, ub)
    if not res.success:
        vf = np.zeros(R); vr = np.zeros(R)
        ef = np.zeros(R); er = np.zeros(R)
    else:
        x = res.x
        vf = x[VF]; vr = x[VR]; ef = x[EF]; er = x[ER]

    return EcFBAResult(reactions=rxns, metabolites=mets, vf=vf, vr=vr, ef=ef, er=er)

def write_flux_table(result: EcFBAResult, rxn_df: pd.DataFrame, out_csv: str):
    flux = pd.DataFrame({
        "reaction_id": result.reactions,
        "vf": result.vf,
        "vr": result.vr,
        "flux_abs": np.abs(result.vf - result.vr),
        "ef": result.ef,
        "er": result.er,
        "enzyme_abs": result.ef + result.er
    })
    merged = rxn_df.merge(flux, on="reaction_id", how="left")
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_csv, index=False)
