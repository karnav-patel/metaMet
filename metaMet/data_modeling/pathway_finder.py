# metaMet/modeling/pathway_finder.py
from __future__ import annotations
import json
import math
from typing import Dict, List, Tuple
from pathlib import Path

import networkx as nx
import pandas as pd

import data_preprocessing.config.config as config
from hybrid.utils_io import safe_read_csv

def build_bipartite_graph(rxn_df: pd.DataFrame, weights_df: pd.DataFrame) -> nx.DiGraph:
    """
    Build a bipartite graph: metabolite -> reaction -> metabolite.
    Edge weights use inverse of reaction weight to prefer high-likelihood steps.
    """
    wmap = dict(zip(weights_df["reaction_id"], weights_df["weight"]))
    G = nx.DiGraph()
    for _, r in rxn_df.iterrows():
        rid = r["reaction_id"]
        w = float(wmap.get(rid, 0.05))
        cost = 1.0 / max(w, 1e-6)  # smaller cost for higher weight
        for m in r["educts"]:
            G.add_edge(f"M::{m}", f"R::{rid}", weight=cost)
        for m in r["products"]:
            G.add_edge(f"R::{rid}", f"M::{m}", weight=cost)
    return G

def shortest_path_between(
    rxn_df: pd.DataFrame,
    weights_df: pd.DataFrame,
    sources: List[str],
    targets: List[str],
    out_json: str
):
    """
    Find shortest paths (Dijkstra) w.r.t cost=1/weight through metabolite-reaction bipartite graph.
    sources/targets are metabolite IDs exactly as used in reaction_index (CUR::ATP, TXT::ethanol, SMI::...).
    """
    G = build_bipartite_graph(rxn_df, weights_df)
    paths = []
    for s in sources:
        sN = f"M::{s}"
        for t in targets:
            tN = f"M::{t}"
            if sN in G and tN in G:
                try:
                    p = nx.shortest_path(G, sN, tN, weight="weight")
                    paths.append({"source": s, "target": t, "path": p})
                except nx.NetworkXNoPath:
                    pass
    Path(out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(out_json).write_text(json.dumps(paths, indent=2))
    return paths
