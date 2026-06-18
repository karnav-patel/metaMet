# metaMet/run_hybrid.py
from __future__ import annotations
import json
from pathlib import Path

import data_preprocessing.config.config as config
from hybrid.aggregate_kcat import build_kcat_aggregate
from hybrid.build_reaction_index import build_reaction_index
from hybrid.reaction_weights import make_reaction_weights
from ec_fba import solve_ecfba, write_flux_table
from pathway_finder import shortest_path_between
from hybrid.utils_io import safe_read_csv

def main(
    seeds=None,
    targets=None,
    enzyme_budget=None
):
    print("Step 1/5: Aggregating kcat (measured + predicted + mapping)...")
    kcat_overview = build_kcat_aggregate()
    print(f"  -> wrote: {config.kcat_aggregate_csv}")
    print(f"  -> wrote: {config.merged_overview}")
    print(f"  -> ECs with kcat: {config.merged_ec_numbers_with_kcat}")
    print(f"  -> ECs missing kcat: {config.predictkcat_missing_kcat_ids_output_path}")

    print("Step 2/5: Building unified reaction index...")
    rxn_df = build_reaction_index()
    print(f"  -> wrote: {config.reaction_index_csv} (n={len(rxn_df)})")

    print("Step 3/5: Creating reaction weights from kcat evidence...")
    weights_df = make_reaction_weights(config.kcat_aggregate_csv, config.reaction_index_csv)
    print(f"  -> wrote: {config.hybrid_weights_csv}")

    print("Step 4/5: Solving enzyme-constrained FBA (hybrid objective)...")
    seeds = seeds if seeds is not None else config.DEFAULT_SEEDS
    targets = targets if targets is not None else config.DEFAULT_TARGETS
    E = float(enzyme_budget) if enzyme_budget is not None else float(config.E_TOTAL)

    result = solve_ecfba(rxn_df, weights_df, seeds=seeds, targets=targets, E_total=E)
    write_flux_table(result, rxn_df, config.ecfba_fluxes_csv)
    print(f"  -> wrote: {config.ecfba_fluxes_csv}")

    print("Step 5/5: Pathway suggestions (shortest high-likelihood routes, optional)...")
    # You can specify explicit metabolite IDs. If none, we skip unless both lists provided.
    if seeds and targets:
        paths = shortest_path_between(rxn_df, weights_df, seeds, targets, config.active_pathways_json)
        print(f"  -> wrote: {config.active_pathways_json} (found {len(paths)} paths)")
    else:
        print("  -> Skipped (provide seeds & targets to activate).")

    print("Done.")

if __name__ == "__main__":
    # Example: edit defaults here or pass via CLI (customize as needed)
    main(
        seeds=["CUR::NAD+", "TXT::primary alcohol"],
        targets=["TXT::aldehyde", "CUR::NADH"],
        enzyme_budget=1.0,
    )
