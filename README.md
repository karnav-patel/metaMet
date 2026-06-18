# ModelYML

Pipeline workspace for building metaMet + GECKO enzyme-constrained models from a protein FASTA.

## Repository layout

- `run_pipeline.py` — the only root-level driver you need to run.
- `config/` — TOML runtime configuration files.
- `inputs/` — user-provided FASTA inputs.
- `src/modelyml/` — Python pipeline source code.
- `matlab/` — MATLAB adapters, helpers, and pipeline scripts.
- `artifacts/` — generated models, adapters, reports, logs, and collected snapshots.
- `archive/oldfiles/` — archived reference outputs used for comparisons.
- `docs/` — workflow notes and static documentation.
- `external/` — auto-managed third-party dependencies.
- `metaMet/` — local metaMet data and modeling assets.

## GECKO handling

The pipeline no longer depends on a hard-coded local GECKO path.

- Default checkout location: `external/GECKO-main`
- If the checkout is missing, the pipeline clones it automatically.
- You can override the location in `config/pipeline.toml`.
- You can override the source repo and ref in `config/pipeline.toml`.

## Main entry point

- `run_pipeline.py` reads `config/pipeline.toml` and drives light, full, or all variants.

## Typical usage

```bash
python run_pipeline.py --config config/pipeline.toml
```

With the default config, that command:
- runs light + full variants
- runs both gapfilled and non-gapfilled cases
- refreshes `artifacts/collections/`
- keeps generated outputs out of the repository root

## Notes

- MATLAB scripts resolve the workspace root, inputs directory, artifacts directory, and GECKO location from environment variables set by the Python runners.
- Generated drafts, models, logs, adapter caches, and collection snapshots are written under `artifacts/`.
- Generated models, logs, cloned GECKO sources, and collected outputs are excluded by `.gitignore`.
