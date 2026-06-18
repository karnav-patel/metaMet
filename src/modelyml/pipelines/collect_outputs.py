from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from modelyml.pipeline_runtime import (  # type: ignore[no-redef]
    DEFAULT_GECKO_DIR,
    DEFAULT_GECKO_REPO,
    gecko_tutorial_output_dir,
    get_workspace_layout,
    infer_workspace_root,
)


INPUT_FILES = [
    Path("inputs/genome.faa"),
    Path("metaMet/data/processed/overview/kcat_aggregate.csv"),
    Path("metaMet/data/raw/mapping_ec_numbers/mapping_ec_number_old_new.csv"),
]


class PipelineError(RuntimeError):
    pass


def parse_args() -> argparse.Namespace:
    root = infer_workspace_root(Path(__file__))
    layout = get_workspace_layout(root)
    parser = argparse.ArgumentParser(
        description=(
            "Run light/full GECKO pipelines with and without CarveMe gapfilling, then collect only the "
            "runtime input files and generated outputs into artifacts/collections."
        )
    )
    parser.add_argument("--workspace", type=Path, default=root, help="Workspace root.")
    parser.add_argument("--python-cmd", default=sys.executable, help="Python executable used to launch the pipelines.")
    parser.add_argument("--gecko-dir", type=Path, default=DEFAULT_GECKO_DIR, help="Where the GECKO checkout should live. Defaults to external/GECKO-main inside the workspace.")
    parser.add_argument("--gecko-repo", default=DEFAULT_GECKO_REPO, help="Git URL used when GECKO needs to be cloned.")
    parser.add_argument("--gecko-ref", default=None, help="Optional git branch or tag to clone for GECKO.")
    parser.add_argument("--matlab-cmd", default=None, help="Optional MATLAB executable path forwarded to the pipelines.")
    parser.add_argument("--carve-cmd", default=None, help="Optional CarveMe executable forwarded to the pipelines.")
    parser.add_argument("--matlab-mode", choices=("batch", "r"), default=None, help="Optional MATLAB mode forwarded to the pipelines.")
    parser.add_argument("--gapfill-media", default="M9", help="Gapfilling medium to use for the *with-gapfilling* runs.")
    parser.add_argument("--output-root", type=Path, default=layout.collections_dir, help="Folder where the four collected cases will be created.")
    parser.add_argument("--resume", action="store_true", help="Forward resume mode to the underlying pipelines.")
    parser.add_argument("--dry-run", action="store_true", help="Print the commands and planned copies without executing them.")
    return parser.parse_args()


def ensure_exists(path: Path, what: str) -> None:
    if not path.exists():
        raise PipelineError(f"Missing {what}: {path}")


def run_command(command: list[str], cwd: Path, dry_run: bool) -> None:
    print("$ " + " ".join(shlex_quote(p) for p in command), flush=True)
    if dry_run:
        return
    completed = subprocess.run(command, cwd=cwd)
    if completed.returncode != 0:
        raise PipelineError(f"Command failed with exit code {completed.returncode}: {' '.join(command)}")


def shlex_quote(text: str) -> str:
    if not text:
        return "''"
    if all(ch.isalnum() or ch in "._/-" for ch in text):
        return text
    return "'" + text.replace("'", "'\"'\"'") + "'"


def make_case_dirs(root: Path, case_name: str, dry_run: bool) -> tuple[Path, Path]:
    input_dir = root / case_name / "input"
    output_dir = root / case_name / "output"
    if dry_run:
        print(f"Would create {input_dir}")
        print(f"Would create {output_dir}")
    else:
        input_dir.mkdir(parents=True, exist_ok=True)
        output_dir.mkdir(parents=True, exist_ok=True)
    return input_dir, output_dir


def copy_inputs(workspace: Path, destination: Path, dry_run: bool) -> None:
    for rel in INPUT_FILES:
        src = workspace / rel
        ensure_exists(src, "input file")
        dst = destination / src.name
        print(f"Copy {src} -> {dst}")
        if not dry_run:
            shutil.copy2(src, dst)


def build_light_output_files(workspace: Path, gecko_dir: Path) -> list[Path]:
    layout = get_workspace_layout(workspace)
    return [
        layout.sanitized_fasta,
        layout.sanitized_fasta_map,
        layout.draft_model_xml,
        layout.draft_model_yaml,
        layout.light_adapter_data_dir / "uniprot.tsv",
        layout.rxn_to_ec_csv,
        layout.light_adapter_data_dir / "customKcats.tsv",
        layout.light_ec_model,
        layout.light_ec_model_light,
        layout.light_ec_model_kcat,
        gecko_dir / "tutorials" / "light_ecModel" / "protocol_light.m",
    ]


def build_full_output_files(workspace: Path, gecko_dir: Path) -> list[Path]:
    layout = get_workspace_layout(workspace)
    tutorial_dir = gecko_tutorial_output_dir(gecko_dir)
    return [
        layout.sanitized_fasta,
        layout.sanitized_fasta_map,
        layout.draft_model_xml,
        layout.draft_model_yaml,
        layout.full_adapter_data_dir / "uniprot.tsv",
        layout.full_adapter_data_dir / "uniprotConversion.tsv",
        layout.full_adapter_data_dir / "pseudoRxns.tsv",
        layout.full_adapter_data_dir / "required_inputs_summary.txt",
        layout.full_adapter_data_dir / "customKcats.tsv",
        layout.full_adapter_data_dir / "ComplexPortal.json",
        layout.full_adapter_data_dir / "kegg.tsv",
        layout.rxn_to_ec_csv,
        layout.full_ec_model,
        layout.full_ec_model_kcat,
        gecko_dir / "tutorials" / "full_ecModel" / "protocol_full.m",
        tutorial_dir / "full_model_summary.tsv",
        tutorial_dir / "full_model_flux_comparison.tsv",
        tutorial_dir / "full_model_ec_number_audit.tsv",
        tutorial_dir / "full_model_kcat_distribution.tsv",
        tutorial_dir / "full_model_top_enzyme_usage.tsv",
        tutorial_dir / "full_model_overview.pdf",
    ]


def copy_outputs(workspace: Path, files: list[Path], destination: Path, log_path: Path, dry_run: bool) -> None:
    for rel in files:
        src = rel if rel.is_absolute() else (workspace / rel)
        if not src.exists():
            continue
        dst = destination / src.name
        print(f"Copy {src} -> {dst}")
        if not dry_run:
            shutil.copy2(src, dst)
    dst_log = destination / log_path.name
    if log_path.resolve() == dst_log.resolve():
        print(f"Log already in place: {dst_log}")
        return
    print(f"Copy {log_path} -> {dst_log}")
    if not dry_run:
        shutil.copy2(log_path, dst_log)


def build_pipeline_command(
    workspace: Path,
    python_cmd: str,
    script_name: str,
    log_path: Path,
    gapfill_media: str | None,
    matlab_cmd: str | None,
    carve_cmd: str | None,
    matlab_mode: str | None,
    gecko_dir: Path,
    gecko_repo: str,
    gecko_ref: str | None,
    resume: bool,
) -> list[str]:
    command = [python_cmd, str(script_name), "--workspace", str(workspace), "--log-file", str(log_path)]
    command.extend(["--gecko-dir", str(gecko_dir), "--gecko-repo", gecko_repo])
    if gecko_ref:
        command.extend(["--gecko-ref", gecko_ref])
    if gapfill_media:
        command.extend(["--gapfill-media", gapfill_media])
    if matlab_cmd:
        command.extend(["--matlab-cmd", matlab_cmd])
    if carve_cmd:
        command.extend(["--carve-cmd", carve_cmd])
    if matlab_mode:
        command.extend(["--matlab-mode", matlab_mode])
    if resume:
        command.append("--resume")
    return command


def run_case(
    workspace: Path,
    output_root: Path,
    case_name: str,
    script_name: str,
    output_files: list[Path],
    gapfill_media: str | None,
    python_cmd: str,
    matlab_cmd: str | None,
    carve_cmd: str | None,
    matlab_mode: str | None,
    gecko_dir: Path,
    gecko_repo: str,
    gecko_ref: str | None,
    resume: bool,
    dry_run: bool,
) -> None:
    print(f"\n=== {case_name} ===", flush=True)
    input_dir, output_dir = make_case_dirs(output_root, case_name, dry_run)
    copy_inputs(workspace, input_dir, dry_run)

    log_path = output_dir / ("pipeline_run.log")
    command = build_pipeline_command(
        workspace=workspace,
        python_cmd=python_cmd,
        script_name=script_name,
        log_path=log_path,
        gapfill_media=gapfill_media,
        matlab_cmd=matlab_cmd,
        carve_cmd=carve_cmd,
        matlab_mode=matlab_mode,
        gecko_dir=gecko_dir,
        gecko_repo=gecko_repo,
        gecko_ref=gecko_ref,
        resume=resume,
    )
    run_command(command, workspace, dry_run)
    copy_outputs(workspace, output_files, output_dir, log_path, dry_run)


def main() -> int:
    args = parse_args()
    workspace = args.workspace.resolve()
    output_root = args.output_root.resolve() if args.output_root.is_absolute() else (workspace / args.output_root).resolve()

    ensure_exists(workspace, "workspace")
    layout = get_workspace_layout(workspace)
    ensure_exists(layout.light_pipeline_script, "light pipeline script")
    ensure_exists(layout.full_pipeline_script, "full pipeline script")

    light_outputs = build_light_output_files(workspace, (workspace / args.gecko_dir).resolve() if not args.gecko_dir.is_absolute() else args.gecko_dir.resolve())
    full_outputs = build_full_output_files(workspace, (workspace / args.gecko_dir).resolve() if not args.gecko_dir.is_absolute() else args.gecko_dir.resolve())

    cases = [
        ("Light-output-no-gapfilling", layout.light_pipeline_script, light_outputs, None),
        ("Light-output-with-gapfilling", layout.light_pipeline_script, light_outputs, args.gapfill_media),
        ("full-output-no-gapfilling", layout.full_pipeline_script, full_outputs, None),
        ("full-output-with-gapfilling", layout.full_pipeline_script, full_outputs, args.gapfill_media),
    ]

    for case_name, script_name, output_files, gapfill_media in cases:
        run_case(
            workspace=workspace,
            output_root=output_root,
            case_name=case_name,
            script_name=script_name,
            output_files=output_files,
            gapfill_media=gapfill_media,
            python_cmd=args.python_cmd,
            matlab_cmd=args.matlab_cmd,
            carve_cmd=args.carve_cmd,
            matlab_mode=args.matlab_mode,
            gecko_dir=args.gecko_dir,
            gecko_repo=args.gecko_repo,
            gecko_ref=args.gecko_ref,
            resume=args.resume,
            dry_run=args.dry_run,
        )

    print(f"\nCollected outputs under {output_root}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except PipelineError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
