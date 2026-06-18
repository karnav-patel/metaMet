from __future__ import annotations

import argparse
import os
import selectors
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from time import monotonic
from typing import Iterable, Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from modelyml.pipeline_runtime import (  # type: ignore[no-redef]
    DEFAULT_GECKO_DIR,
    DEFAULT_GECKO_REPO,
    RuntimeSetupError,
    build_runtime_env,
    ensure_gecko_checkout,
    ensure_runtime_layout,
    get_workspace_layout,
    infer_workspace_root,
    normalize_relpath,
    resolve_gecko_dir,
)


@dataclass(frozen=True)
class Step:
    number: int
    name: str
    kind: str
    command: Sequence[str] | None = None
    matlab_script: str | None = None
    expected_outputs: tuple[str, ...] = ()


class PipelineError(RuntimeError):
    pass


SANITIZED_FASTA = "genome_for_carveme.faa"
SANITIZED_FASTA_MAP = "genome_for_carveme_header_map.tsv"
NON_GAPFILL_MODEL = "model_nogapfill.xml"
DEFAULT_GAPFILL_REPORT = "gapfill_compare_report.md"


def parse_args() -> argparse.Namespace:
    root = infer_workspace_root(Path(__file__))
    parser = argparse.ArgumentParser(
        description=(
            "Run the metaMet ↔ GECKO ecModel pipeline from genome.faa to ecModel_kcat.yml. "
            "The script orchestrates CarveMe plus the MATLAB scripts already present in this workspace."
        )
    )
    parser.add_argument(
        "--workspace",
        type=Path,
        default=root,
        help="Workspace root containing genome.faa, the MATLAB scripts, and metaMet.",
    )
    parser.add_argument(
        "--gecko-dir",
        type=Path,
        default=DEFAULT_GECKO_DIR,
        help="Where the GECKO checkout should live. Defaults to external/GECKO-main inside the workspace.",
    )
    parser.add_argument(
        "--gecko-repo",
        default=DEFAULT_GECKO_REPO,
        help="Git URL used when GECKO needs to be cloned.",
    )
    parser.add_argument(
        "--gecko-ref",
        default=None,
        help="Optional git branch or tag to clone for GECKO.",
    )
    parser.add_argument(
        "--matlab-cmd",
        default=os.environ.get("MATLAB_CMD"),
        help="Path to the MATLAB executable. Defaults to $MATLAB_CMD or auto-detection.",
    )
    parser.add_argument(
        "--carve-cmd",
        default=os.environ.get("CARVE_CMD") or "carve",
        help="CarveMe executable or command name. Defaults to 'carve'.",
    )
    parser.add_argument(
        "--gapfill-media",
        default=None,
        help=(
            "Optional comma-separated CarveMe media list for gap-filling during reconstruction "
            "(for example: M9 or M9,LB)."
        ),
    )
    parser.add_argument(
        "--init-media",
        default=None,
        help=(
            "Optional CarveMe medium used to initialize exchange bounds in the draft model "
            "without changing gap-filling behavior unless --gapfill-media is also set."
        ),
    )
    parser.add_argument(
        "--compare-gapfill",
        action="store_true",
        help=(
            "When used with --gapfill-media, also build a non-gapfilled CarveMe draft and write "
            "a Markdown comparison report."
        ),
    )
    parser.add_argument(
        "--gapfill-report",
        type=Path,
        default=Path(DEFAULT_GAPFILL_REPORT),
        help=(
            "Markdown report path for the gapfilled-versus-non-gapfilled comparison. "
            f"Defaults to {DEFAULT_GAPFILL_REPORT} under the workspace."
        ),
    )
    parser.add_argument(
        "--matlab-mode",
        choices=("batch", "r"),
        default="batch",
        help="How to invoke MATLAB. Use 'batch' for modern MATLAB or 'r' for older releases.",
    )
    parser.add_argument(
        "--from-step",
        type=int,
        default=1,
        help="First pipeline step to run (1-6).",
    )
    parser.add_argument(
        "--to-step",
        type=int,
        default=6,
        help="Last pipeline step to run (1-6).",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip a step when all of its expected output files already exist.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the commands that would run without executing them.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=None,
        help="Optional log file path. Defaults to pipeline_logs/<timestamp>.log under the workspace.",
    )
    return parser.parse_args()


def build_steps(
    workspace: Path,
    carve_cmd: str,
    python_cmd: str,
    gapfill_media: str | None = None,
    init_media: str | None = None,
    compare_gapfill: bool = False,
    gapfill_report: Path | None = None,
) -> list[Step]:
    layout = get_workspace_layout(workspace)

    def rel(path: Path) -> str:
        return normalize_relpath(path.relative_to(workspace))

    def build_carve_command(output_name: Path, gapfill: str | None) -> tuple[str, ...]:
        command: list[str] = [carve_cmd, str(layout.sanitized_fasta), "-o", str(output_name)]
        if gapfill:
            command.extend(["--gapfill", gapfill])
        if init_media:
            command.extend(["--init", init_media])
        return tuple(command)

    steps = [
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
    ]

    if compare_gapfill and gapfill_media:
        steps.append(
            Step(
                number=1,
                name="Build non-gapfilled draft GEM with CarveMe",
                kind="command",
                command=build_carve_command(layout.non_gapfill_model_xml, None),
                expected_outputs=(rel(layout.non_gapfill_model_xml),),
            )
        )

    steps.extend(
        [
            Step(
            number=1,
            name="Build draft GEM with CarveMe",
            kind="command",
            command=build_carve_command(layout.draft_model_xml, gapfill_media),
            expected_outputs=(rel(layout.draft_model_xml),),
        ),
        ]
    )

    if compare_gapfill and gapfill_media and gapfill_report is not None:
        steps.append(
            Step(
                number=1,
                name="Write gapfill comparison report",
                kind="command",
                command=(
                    python_cmd,
                    str(layout.compare_analysis_script),
                    "--current-file",
                    str(layout.draft_model_xml),
                    "--reference-file",
                    str(layout.non_gapfill_model_xml),
                    "--kind",
                    "xml",
                    "--report",
                    str(gapfill_report),
                    "--report-title",
                    "Gapfill comparison report",
                    "--report-description",
                    "This report compares the gapfilled CarveMe draft against a non-gapfilled draft built from the same sanitized FASTA input.",
                ),
            )
        )

    steps.extend([
        Step(
            number=2,
            name="Convert SBML to GECKO/RAVEN YAML",
            kind="matlab",
            matlab_script="sbml_to_yaml_roundtrip.m",
            expected_outputs=(rel(layout.draft_model_yaml),),
        ),
        Step(
            number=3,
            name="Prepare adapter and build ecModel scaffold",
            kind="matlab",
            matlab_script="prepare_gecko_adapter.m",
            expected_outputs=(rel(layout.light_adapter_data_dir / "uniprot.tsv"),),
        ),
        Step(
            number=3,
            name="Build ecModel scaffold",
            kind="matlab",
            matlab_script="build_ecmodel_gecko.m",
            expected_outputs=(rel(layout.light_ec_model),),
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
            name="Apply metaMet kcats to ecModel",
            kind="matlab",
            matlab_script="apply_kcat_from_metMet.m",
            expected_outputs=(
                rel(layout.light_adapter_data_dir / "customKcats.tsv"),
                rel(layout.light_ec_model_kcat),
                "metaMet/data_modeling/models/ecModel_kcat.yml",
            ),
        ),
        Step(
            number=6,
            name="Run QC checks",
            kind="matlab",
            matlab_script="qc_gecko_models.m",
            expected_outputs=(rel(layout.light_ec_model_kcat),),
        ),
    ])

    return steps


def ensure_workspace_layout(workspace: Path, selected_steps: Iterable[Step]) -> None:
    layout = get_workspace_layout(workspace)
    required = [
        layout.genome_faa,
        layout.sanitize_fasta_script,
        "metaMet",
        layout.matlab_adapters_dir / "CarveMeModelAdapter.m",
        layout.matlab_scripts_dir / "prepare_gecko_adapter.m",
        layout.matlab_scripts_dir / "sbml_to_yaml_roundtrip.m",
        layout.matlab_scripts_dir / "build_ecmodel_gecko.m",
        layout.matlab_scripts_dir / "export_rxn_to_ec.m",
        layout.matlab_scripts_dir / "apply_kcat_from_metMet.m",
        layout.matlab_scripts_dir / "qc_gecko_models.m",
    ]
    if any(step.command and any("compare_model_analysis.py" in part for part in step.command) for step in selected_steps):
        required.append(layout.compare_analysis_script)
    missing = [
        str(item if isinstance(item, Path) else (workspace / item))
        for item in required
        if not (item if isinstance(item, Path) else (workspace / item)).exists()
    ]
    if missing:
        joined = ", ".join(missing)
        raise PipelineError(f"Workspace is missing required pipeline files/directories: {joined}")

    if any(step.number == 5 for step in selected_steps):
        kcat_table = workspace / "metaMet/data/processed/overview/kcat_aggregate.csv"
        if not kcat_table.is_file():
            raise PipelineError(f"Required kcat table not found: {kcat_table}")


def default_log_file(workspace: Path) -> Path:
    log_dir = get_workspace_layout(workspace).logs_dir
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return log_dir / f"metaMet_gecko_pipeline_{stamp}.log"


def quote_for_matlab(path: Path | str) -> str:
    return str(path).replace("\\", "/").replace("'", "''")


def resolve_matlab_command(matlab_cmd: str | None) -> str:
    candidates: list[str] = []
    if matlab_cmd:
        candidates.append(matlab_cmd)
    env_cmd = os.environ.get("MATLAB_CMD")
    if env_cmd and env_cmd not in candidates:
        candidates.append(env_cmd)
    which_matlab = shutil.which("matlab")
    if which_matlab:
        candidates.append(which_matlab)

    applications = Path("/Applications")
    if applications.is_dir():
        app_bins = sorted(applications.glob("MATLAB*.app/bin/matlab"), reverse=True)
        candidates.extend(str(path) for path in app_bins)

    for candidate in candidates:
        if not candidate:
            continue
        resolved = shutil.which(candidate) if os.sep not in candidate else candidate
        if resolved and Path(resolved).exists():
            return resolved

    raise PipelineError(
        "MATLAB executable was not found. Set --matlab-cmd or the MATLAB_CMD environment variable."
    )


def resolve_command_path(command: str) -> str:
    if os.sep in command:
        if Path(command).exists():
            return command
        raise PipelineError(f"Command path does not exist: {command}")

    resolved = shutil.which(command)
    if resolved:
        return resolved
    raise PipelineError(
        f"Command '{command}' was not found on PATH. Install it or pass an explicit path."
    )


def step_outputs_exist(workspace: Path, step: Step) -> bool:
    return bool(step.expected_outputs) and all((workspace / rel).exists() for rel in step.expected_outputs)


def emit(message: str, log_handle) -> None:
    print(message, flush=True)
    log_handle.write(message + "\n")
    log_handle.flush()


def format_command(command: Sequence[str]) -> str:
    return " ".join(shlex.quote(part) for part in command)


def run_process(command: Sequence[str], cwd: Path, log_handle, dry_run: bool, env: dict[str, str] | None = None) -> None:
    emit(f"$ {format_command(command)}", log_handle)
    if dry_run:
        return

    process = subprocess.Popen(
        list(command),
        cwd=cwd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert process.stdout is not None
    selector = selectors.DefaultSelector()
    selector.register(process.stdout, selectors.EVENT_READ)
    exit_grace_started: float | None = None
    open_pipe_notice_emitted = False

    try:
        while True:
            events = selector.select(timeout=0.5)
            if events:
                for key, _ in events:
                    line = key.fileobj.readline()
                    if line == "":
                        try:
                            selector.unregister(key.fileobj)
                        except Exception:
                            pass
                        continue
                    text = line.rstrip("\n")
                    print(text, flush=True)
                    log_handle.write(text + "\n")
                log_handle.flush()
                exit_grace_started = None

            if process.poll() is None:
                continue

            if not selector.get_map():
                break

            if exit_grace_started is None:
                exit_grace_started = monotonic()
                continue

            if monotonic() - exit_grace_started < 1.0:
                continue

            if not open_pipe_notice_emitted:
                emit(
                    "NOTE: Child process exited but its stdout pipe remained open via descendant helper processes; continuing after a short grace period.",
                    log_handle,
                )
                open_pipe_notice_emitted = True
            break
    finally:
        selector.close()
        process.stdout.close()

    process.wait()
    log_handle.flush()
    if process.returncode != 0:
        raise PipelineError(f"Command failed with exit code {process.returncode}: {format_command(command)}")


def build_matlab_command(matlab_exe: str, mode: str, workspace: Path, script_name: str) -> list[str]:
    layout = get_workspace_layout(workspace)
    workspace_q = quote_for_matlab(workspace)
    matlab_root_q = quote_for_matlab(layout.matlab_root)
    script_q = quote_for_matlab(layout.matlab_scripts_dir / script_name)
    body = (
        f"addpath(genpath('{matlab_root_q}')); cd('{workspace_q}'); "
        f"try, run('{script_q}'); "
        f"catch ME, disp(getReport(ME,'extended')); exit(1); end"
    )
    if mode == "batch":
        return [matlab_exe, "-batch", body]
    return [matlab_exe, "-nodesktop", "-nosplash", "-r", body + "; exit(0);"]


def run_step(step: Step, workspace: Path, matlab_exe: str, matlab_mode: str, log_handle, dry_run: bool, env: dict[str, str]) -> None:
    emit(f"\n=== Step {step.number}: {step.name} ===", log_handle)
    if step.kind == "command":
        assert step.command is not None
        command = list(step.command)
        command[0] = resolve_command_path(command[0])
        run_process(command, workspace, log_handle, dry_run=dry_run, env=env)
        return

    if step.kind == "matlab":
        assert step.matlab_script is not None
        script_path = get_workspace_layout(workspace).matlab_scripts_dir / step.matlab_script
        if not script_path.is_file():
            raise PipelineError(f"MATLAB script not found: {script_path}")
        command = build_matlab_command(matlab_exe, matlab_mode, workspace, step.matlab_script)
        run_process(command, workspace, log_handle, dry_run=dry_run, env=env)
        return

    raise PipelineError(f"Unsupported step type: {step.kind}")


def verify_outputs(workspace: Path, selected_steps: Iterable[Step], dry_run: bool) -> None:
    if dry_run:
        return

    missing: list[str] = []
    for step in selected_steps:
        for rel in step.expected_outputs:
            path = workspace / rel
            if not path.exists():
                missing.append(str(path))

    if missing:
        raise PipelineError("Expected output files were not produced:\n- " + "\n- ".join(missing))

    final_output = get_workspace_layout(workspace).light_ec_model_kcat
    if any(step.number >= 5 for step in selected_steps) and not final_output.is_file():
        raise PipelineError(f"Final deliverable was not produced: {final_output}")


def select_steps(all_steps: Sequence[Step], from_step: int, to_step: int) -> list[Step]:
    max_step = max((step.number for step in all_steps), default=0)
    if from_step < 1 or to_step > max_step or from_step > to_step:
        raise PipelineError(
            f"Step range must satisfy 1 <= --from-step <= --to-step <= {max_step}."
        )
    return [step for step in all_steps if from_step <= step.number <= to_step]


def maybe_skip_step(step: Step, workspace: Path, resume: bool, log_handle) -> bool:
    if resume and step_outputs_exist(workspace, step):
        emit(f"Skipping step {step.number} ({step.name}) because outputs already exist.", log_handle)
        return True
    return False


def main() -> int:
    args = parse_args()
    workspace = args.workspace.resolve()
    layout = ensure_runtime_layout(workspace, dry_run=args.dry_run)
    gecko_dir = resolve_gecko_dir(workspace, args.gecko_dir)
    log_path = args.log_file.resolve() if args.log_file else default_log_file(workspace)
    gapfill_report = (layout.reports_dir / args.gapfill_report).resolve() if not args.gapfill_report.is_absolute() else args.gapfill_report.resolve()

    if not workspace.is_dir():
        raise PipelineError(f"Workspace directory does not exist: {workspace}")
    if args.compare_gapfill and not args.gapfill_media:
        raise PipelineError("--compare-gapfill requires --gapfill-media.")

    log_path.parent.mkdir(parents=True, exist_ok=True)
    gapfill_report.parent.mkdir(parents=True, exist_ok=True)

    steps = select_steps(
        build_steps(
            workspace,
            args.carve_cmd,
            sys.executable,
            gapfill_media=args.gapfill_media,
            init_media=args.init_media,
            compare_gapfill=args.compare_gapfill,
            gapfill_report=gapfill_report,
        ),
        args.from_step,
        args.to_step,
    )
    ensure_workspace_layout(workspace, steps)
    matlab_exe = ""
    runtime_env = build_runtime_env(workspace, gecko_dir)

    with log_path.open("w", encoding="utf-8") as log_handle:
        emit("metaMet ↔ GECKO pipeline runner", log_handle)
        emit(f"Workspace: {workspace}", log_handle)
        emit(f"GECKO checkout: {gecko_dir}", log_handle)
        emit(f"Log file: {log_path}", log_handle)
        emit(f"Dry run: {'yes' if args.dry_run else 'no'}", log_handle)
        emit(f"Resume mode: {'yes' if args.resume else 'no'}", log_handle)
        if args.compare_gapfill and args.gapfill_media:
            emit(f"Gapfill comparison report: {gapfill_report}", log_handle)

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
        emit("\nPipeline completed successfully.", log_handle)
        if any(step.number >= 5 for step in steps):
            emit(f"Final deliverable: {layout.light_ec_model_kcat}", log_handle)
        if args.compare_gapfill and args.gapfill_media and not args.dry_run:
            emit(f"Gapfill draft-model report: {gapfill_report}", log_handle)

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
