from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path
from typing import Iterable, Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from modelyml.pipelines.light_pipeline import (  # type: ignore[no-redef]
    PipelineError,
    SANITIZED_FASTA,
    SANITIZED_FASTA_MAP,
    Step,
    default_log_file,
    emit,
    maybe_skip_step,
    resolve_matlab_command,
    run_step,
    select_steps,
)
from modelyml.pipeline_runtime import (  # type: ignore[no-redef]
    DEFAULT_GECKO_DIR,
    DEFAULT_GECKO_REPO,
    RuntimeSetupError,
    build_runtime_env,
    ensure_gecko_checkout,
    ensure_runtime_layout,
    gecko_tutorial_output_dir,
    get_workspace_layout,
    infer_workspace_root,
    normalize_relpath,
    resolve_gecko_dir,
)


FULL_TUTORIAL_OUTPUT = DEFAULT_GECKO_DIR / "tutorials" / "full_ecModel" / "output"


def parse_args() -> argparse.Namespace:
    root = infer_workspace_root(Path(__file__))
    parser = argparse.ArgumentParser(
        description=(
            "Run the metaMet ↔ GECKO full-model pipeline from genome.faa to ecModel_full_kcat.yml. "
            "This pipeline prepares the extra adapter inputs required by the full GECKO formulation."
        )
    )
    parser.add_argument("--workspace", type=Path, default=root, help="Workspace root.")
    parser.add_argument("--gecko-dir", type=Path, default=DEFAULT_GECKO_DIR, help="Where the GECKO checkout should live. Defaults to external/GECKO-main inside the workspace.")
    parser.add_argument("--gecko-repo", default=DEFAULT_GECKO_REPO, help="Git URL used when GECKO needs to be cloned.")
    parser.add_argument("--gecko-ref", default=None, help="Optional git branch or tag to clone for GECKO.")
    parser.add_argument("--matlab-cmd", default=os.environ.get("MATLAB_CMD"), help="Path to MATLAB executable.")
    parser.add_argument("--carve-cmd", default=os.environ.get("CARVE_CMD") or "carve", help="CarveMe executable or command name.")
    parser.add_argument("--gapfill-media", default=None, help="Optional comma-separated CarveMe media list for gap-filling during reconstruction.")
    parser.add_argument("--init-media", default=None, help="Optional CarveMe medium used to initialize exchange bounds in the draft model.")
    parser.add_argument("--matlab-mode", choices=("batch", "r"), default="batch", help="How to invoke MATLAB.")
    parser.add_argument("--from-step", type=int, default=1, help="First pipeline step to run (1-7).")
    parser.add_argument("--to-step", type=int, default=7, help="Last pipeline step to run (1-7).")
    parser.add_argument("--resume", action="store_true", help="Skip a step when all expected outputs already exist.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing them.")
    parser.add_argument("--log-file", type=Path, default=None, help="Optional log file path.")
    return parser.parse_args()


def build_steps(workspace: Path, carve_cmd: str, python_cmd: str, gapfill_media: str | None = None, init_media: str | None = None) -> list[Step]:
    layout = get_workspace_layout(workspace)

    def rel(path: Path) -> str:
        return normalize_relpath(path.relative_to(workspace))

    carve_command: list[str] = [carve_cmd, str(layout.sanitized_fasta), "-o", str(layout.draft_model_xml)]
    if gapfill_media:
        carve_command.extend(["--gapfill", gapfill_media])
    if init_media:
        carve_command.extend(["--init", init_media])

    return [
        Step(
            number=1,
            name="Sanitize FASTA headers for CarveMe",
            kind="command",
            command=(
                python_cmd,
                str(layout.sanitize_fasta_script),
                "--input",
                str(layout.genome_faa),
                "--output",
                str(layout.sanitized_fasta),
                "--map-out",
                str(layout.sanitized_fasta_map),
            ),
            expected_outputs=(rel(layout.sanitized_fasta), rel(layout.sanitized_fasta_map)),
        ),
        Step(
            number=1,
            name="Build draft GEM with CarveMe",
            kind="command",
            command=tuple(carve_command),
            expected_outputs=(rel(layout.draft_model_xml),),
        ),
        Step(
            number=2,
            name="Convert SBML to GECKO/RAVEN YAML",
            kind="matlab",
            matlab_script="sbml_to_yaml_roundtrip.m",
            expected_outputs=(rel(layout.draft_model_yaml),),
        ),
        Step(
            number=3,
            name="Prepare full-model GECKO adapter inputs",
            kind="matlab",
            matlab_script="prepare_gecko_full_adapter.m",
            expected_outputs=(
                rel(layout.full_adapter_data_dir / "uniprot.tsv"),
                rel(layout.full_adapter_data_dir / "uniprotConversion.tsv"),
                rel(layout.full_adapter_data_dir / "pseudoRxns.tsv"),
            ),
        ),
        Step(
            number=3,
            name="Build full GECKO ecModel scaffold",
            kind="matlab",
            matlab_script="build_ecmodel_gecko_full.m",
            expected_outputs=(rel(layout.full_ec_model),),
        ),
        Step(
            number=4,
            name="Export reaction to EC mapping",
            kind="matlab",
            matlab_script="export_rxn_to_ec.m",
            expected_outputs=(rel(layout.rxn_to_ec_csv),),
        ),
        Step(
            number=5,
            name="Apply metaMet kcats to full GECKO model",
            kind="matlab",
            matlab_script="apply_kcat_from_metMet_full.m",
            expected_outputs=(
                rel(layout.full_adapter_data_dir / "customKcats.tsv"),
                rel(layout.full_ec_model_kcat),
                "metaMet/data_modeling/models/ecModel_full_kcat.yml",
            ),
        ),
        Step(
            number=6,
            name="Run QC checks on full GECKO models",
            kind="matlab",
            matlab_script="qc_gecko_full_models.m",
            expected_outputs=(rel(layout.full_ec_model_kcat),),
        ),
        Step(
            number=7,
            name="Generate rich full-model analysis outputs",
            kind="matlab",
            matlab_script="analyze_gecko_full_outputs.m",
            expected_outputs=(
                normalize_relpath(FULL_TUTORIAL_OUTPUT / "full_model_summary.tsv"),
                normalize_relpath(FULL_TUTORIAL_OUTPUT / "full_model_flux_comparison.tsv"),
                normalize_relpath(FULL_TUTORIAL_OUTPUT / "full_model_kcat_distribution.tsv"),
                normalize_relpath(FULL_TUTORIAL_OUTPUT / "full_model_top_enzyme_usage.tsv"),
                normalize_relpath(FULL_TUTORIAL_OUTPUT / "full_model_overview.pdf"),
            ),
        ),
    ]


def ensure_workspace_layout(workspace: Path, selected_steps: Iterable[Step]) -> None:
    layout = get_workspace_layout(workspace)
    required = [
        layout.genome_faa,
        layout.sanitize_fasta_script,
        "metaMet",
        layout.matlab_adapters_dir / "CarveMeFullModelAdapter.m",
        layout.matlab_scripts_dir / "prepare_gecko_full_adapter.m",
        layout.matlab_scripts_dir / "sbml_to_yaml_roundtrip.m",
        layout.matlab_scripts_dir / "build_ecmodel_gecko_full.m",
        layout.matlab_scripts_dir / "export_rxn_to_ec.m",
        layout.matlab_scripts_dir / "apply_kcat_from_metMet_full.m",
        layout.matlab_scripts_dir / "qc_gecko_full_models.m",
        layout.matlab_scripts_dir / "analyze_gecko_full_outputs.m",
    ]
    missing = [
        str(item if isinstance(item, Path) else (workspace / item))
        for item in required
        if not (item if isinstance(item, Path) else (workspace / item)).exists()
    ]
    if missing:
        raise PipelineError(f"Workspace is missing required pipeline files/directories: {', '.join(missing)}")

    if any(step.number == 5 for step in selected_steps):
        kcat_table = workspace / "metaMet/data/processed/overview/kcat_aggregate.csv"
        if not kcat_table.is_file():
            raise PipelineError(f"Required kcat table not found: {kcat_table}")


def verify_outputs(workspace: Path, selected_steps: Sequence[Step], dry_run: bool) -> None:
    if dry_run:
        return

    missing: list[str] = []
    for step in selected_steps:
        for rel in step.expected_outputs:
            if not (workspace / rel).exists():
                missing.append(str(workspace / rel))
    if missing:
        raise PipelineError("Expected output files were not produced:\n- " + "\n- ".join(missing))

    final_output = get_workspace_layout(workspace).full_ec_model_kcat
    if any(step.number >= 5 for step in selected_steps) and not final_output.is_file():
        raise PipelineError(f"Final deliverable was not produced: {final_output}")


def sync_tutorial_output(workspace: Path, gecko_dir: Path, log_path: Path, selected_steps: Sequence[Step], dry_run: bool) -> None:
    if dry_run:
        return

    layout = get_workspace_layout(workspace)
    tutorial_output = gecko_tutorial_output_dir(gecko_dir)
    if not tutorial_output.is_dir():
        return

    files_to_copy = {
        layout.full_ec_model: tutorial_output / "ecModel_full.yml",
        layout.full_ec_model_kcat: tutorial_output / "ecModel_full_kcat.yml",
        layout.rxn_to_ec_csv: tutorial_output / "rxn_to_ec_full.csv",
        layout.full_adapter_data_dir / "customKcats.tsv": tutorial_output / "customKcats_full.tsv",
        log_path: tutorial_output / "full_ecModel_run.log",
    }

    for source, target in files_to_copy.items():
        if source.is_file():
            shutil.copy2(source, target)


def main() -> int:
    args = parse_args()
    workspace = args.workspace.resolve()
    layout = ensure_runtime_layout(workspace, dry_run=args.dry_run)
    gecko_dir = resolve_gecko_dir(workspace, args.gecko_dir)
    log_path = args.log_file.resolve() if args.log_file else default_log_file(workspace)

    if not workspace.is_dir():
        raise PipelineError(f"Workspace directory does not exist: {workspace}")

    log_path.parent.mkdir(parents=True, exist_ok=True)

    steps = select_steps(build_steps(workspace, args.carve_cmd, sys.executable, args.gapfill_media, args.init_media), args.from_step, args.to_step)
    ensure_workspace_layout(workspace, steps)
    matlab_exe = ""
    runtime_env = build_runtime_env(workspace, gecko_dir)

    with log_path.open("w", encoding="utf-8") as log_handle:
        emit("metaMet ↔ GECKO full-model pipeline runner", log_handle)
        emit(f"Workspace: {workspace}", log_handle)
        emit(f"GECKO checkout: {gecko_dir}", log_handle)
        emit(f"Log file: {log_path}", log_handle)
        emit(f"Dry run: {'yes' if args.dry_run else 'no'}", log_handle)
        emit(f"Resume mode: {'yes' if args.resume else 'no'}", log_handle)

        try:
            ensure_gecko_checkout(
                workspace=workspace,
                gecko_dir=gecko_dir,
                repo_url=args.gecko_repo,
                ref=args.gecko_ref,
                logger=lambda message: emit(message, log_handle),
                dry_run=args.dry_run,
            )
        except RuntimeSetupError as exc:
            raise PipelineError(str(exc)) from exc

        for step in steps:
            if maybe_skip_step(step, workspace, args.resume, log_handle):
                continue
            if step.kind == "matlab" and not matlab_exe:
                matlab_exe = resolve_matlab_command(args.matlab_cmd)
                emit(f"MATLAB: {matlab_exe}", log_handle)
            run_step(step, workspace, matlab_exe, args.matlab_mode, log_handle, args.dry_run, runtime_env)

        if not matlab_exe:
            emit("MATLAB: not needed", log_handle)

        verify_outputs(workspace, steps, args.dry_run)
        emit("\nFull-model pipeline completed successfully.", log_handle)
        if any(step.number >= 5 for step in steps):
            emit(f"Final deliverable: {layout.full_ec_model_kcat}", log_handle)

    sync_tutorial_output(workspace, gecko_dir, log_path, steps, args.dry_run)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except PipelineError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
    except RuntimeSetupError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
    except KeyboardInterrupt:
        print("ERROR: Pipeline interrupted by user.", file=sys.stderr)
        raise SystemExit(130)