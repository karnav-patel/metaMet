from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


DEFAULT_GECKO_REPO = "https://github.com/SysBioChalmers/GECKO.git"
DEFAULT_GECKO_DIR = Path("external") / "GECKO-main"
DEFAULT_ARCHIVE_DIR = Path("archive") / "oldfiles"
DEFAULT_INPUTS_DIR = Path("inputs")
DEFAULT_ARTIFACTS_DIR = Path("artifacts")
DEFAULT_LOGS_DIR = DEFAULT_ARTIFACTS_DIR / "logs"
DEFAULT_REPORTS_DIR = DEFAULT_ARTIFACTS_DIR / "reports"
DEFAULT_COLLECTIONS_DIR = DEFAULT_ARTIFACTS_DIR / "collections"
DEFAULT_DRAFTS_DIR = DEFAULT_ARTIFACTS_DIR / "drafts"
DEFAULT_MODELS_DIR = DEFAULT_ARTIFACTS_DIR / "models"
DEFAULT_LIGHT_MODELS_DIR = DEFAULT_MODELS_DIR / "light"
DEFAULT_FULL_MODELS_DIR = DEFAULT_MODELS_DIR / "full"
DEFAULT_ADAPTERS_DIR = DEFAULT_ARTIFACTS_DIR / "adapters"
DEFAULT_LIGHT_ADAPTER_DIR = DEFAULT_ADAPTERS_DIR / "light"
DEFAULT_FULL_ADAPTER_DIR = DEFAULT_ADAPTERS_DIR / "full"


class RuntimeSetupError(RuntimeError):
    pass


Logger = Callable[[str], None] | None


@dataclass(frozen=True)
class WorkspaceLayout:
    workspace: Path

    @property
    def src_dir(self) -> Path:
        return self.workspace / "src"

    @property
    def python_package_dir(self) -> Path:
        return self.src_dir / "modelyml"

    @property
    def python_pipelines_dir(self) -> Path:
        return self.python_package_dir / "pipelines"

    @property
    def python_utils_dir(self) -> Path:
        return self.python_package_dir / "utils"

    @property
    def matlab_root(self) -> Path:
        return self.workspace / "matlab"

    @property
    def matlab_scripts_dir(self) -> Path:
        return self.matlab_root / "scripts"

    @property
    def matlab_helpers_dir(self) -> Path:
        return self.matlab_root / "helpers"

    @property
    def matlab_adapters_dir(self) -> Path:
        return self.matlab_root / "adapters"

    @property
    def inputs_dir(self) -> Path:
        return self.workspace / DEFAULT_INPUTS_DIR

    @property
    def artifacts_dir(self) -> Path:
        return self.workspace / DEFAULT_ARTIFACTS_DIR

    @property
    def logs_dir(self) -> Path:
        return self.workspace / DEFAULT_LOGS_DIR

    @property
    def reports_dir(self) -> Path:
        return self.workspace / DEFAULT_REPORTS_DIR

    @property
    def collections_dir(self) -> Path:
        return self.workspace / DEFAULT_COLLECTIONS_DIR

    @property
    def drafts_dir(self) -> Path:
        return self.workspace / DEFAULT_DRAFTS_DIR

    @property
    def models_dir(self) -> Path:
        return self.workspace / DEFAULT_MODELS_DIR

    @property
    def light_models_dir(self) -> Path:
        return self.workspace / DEFAULT_LIGHT_MODELS_DIR

    @property
    def full_models_dir(self) -> Path:
        return self.workspace / DEFAULT_FULL_MODELS_DIR

    @property
    def adapters_dir(self) -> Path:
        return self.workspace / DEFAULT_ADAPTERS_DIR

    @property
    def light_adapter_dir(self) -> Path:
        return self.workspace / DEFAULT_LIGHT_ADAPTER_DIR

    @property
    def full_adapter_dir(self) -> Path:
        return self.workspace / DEFAULT_FULL_ADAPTER_DIR

    @property
    def archive_dir(self) -> Path:
        return self.workspace / DEFAULT_ARCHIVE_DIR

    @property
    def metaMet_models_dir(self) -> Path:
        return self.workspace / "metaMet" / "data_modeling" / "models"

    @property
    def genome_faa(self) -> Path:
        return self.inputs_dir / "genome.faa"

    @property
    def example_genome_faa(self) -> Path:
        return self.inputs_dir / "genome_ABC.faa"

    @property
    def sanitized_fasta(self) -> Path:
        return self.drafts_dir / "genome_for_carveme.faa"

    @property
    def sanitized_fasta_map(self) -> Path:
        return self.drafts_dir / "genome_for_carveme_header_map.tsv"

    @property
    def draft_model_xml(self) -> Path:
        return self.drafts_dir / "model.xml"

    @property
    def draft_model_yaml(self) -> Path:
        return self.drafts_dir / "model.yml"

    @property
    def non_gapfill_model_xml(self) -> Path:
        return self.drafts_dir / "model_nogapfill.xml"

    @property
    def rxn_to_ec_csv(self) -> Path:
        return self.drafts_dir / "rxn_to_ec.csv"

    @property
    def light_ec_model(self) -> Path:
        return self.light_models_dir / "ecModel.yml"

    @property
    def light_ec_model_light(self) -> Path:
        return self.light_models_dir / "ecModel_light.yml"

    @property
    def light_ec_model_kcat(self) -> Path:
        return self.light_models_dir / "ecModel_kcat.yml"

    @property
    def full_ec_model(self) -> Path:
        return self.full_models_dir / "ecModel_full.yml"

    @property
    def full_ec_model_kcat(self) -> Path:
        return self.full_models_dir / "ecModel_full_kcat.yml"

    @property
    def light_adapter_data_dir(self) -> Path:
        return self.light_adapter_dir / "data"

    @property
    def full_adapter_data_dir(self) -> Path:
        return self.full_adapter_dir / "data"

    @property
    def sanitize_fasta_script(self) -> Path:
        return self.python_utils_dir / "sanitize_fasta_for_carveme.py"

    @property
    def compare_analysis_script(self) -> Path:
        return self.python_utils_dir / "compare_model_analysis.py"

    @property
    def light_pipeline_script(self) -> Path:
        return self.python_pipelines_dir / "light_pipeline.py"

    @property
    def full_pipeline_script(self) -> Path:
        return self.python_pipelines_dir / "full_pipeline.py"


def normalize_relpath(path: Path | str) -> str:
    return str(Path(path).as_posix())


def infer_workspace_root(current_file: Path, levels_up: int = 3) -> Path:
    path = current_file.resolve()
    for _ in range(levels_up):
        path = path.parent
    return path


def get_workspace_layout(workspace: Path) -> WorkspaceLayout:
    return WorkspaceLayout(workspace=workspace.resolve())


def resolve_workspace_path(workspace: Path, path: Path) -> Path:
    return path.resolve() if path.is_absolute() else (workspace / path).resolve()


def resolve_gecko_dir(workspace: Path, gecko_dir: Path | None = None) -> Path:
    candidate = gecko_dir or DEFAULT_GECKO_DIR
    return resolve_workspace_path(workspace, candidate)


def gecko_tutorial_output_dir(gecko_dir: Path) -> Path:
    return gecko_dir / "tutorials" / "full_ecModel" / "output"


def ensure_runtime_layout(workspace: Path, dry_run: bool = False) -> WorkspaceLayout:
    layout = get_workspace_layout(workspace)
    for path in (
        workspace / DEFAULT_GECKO_DIR.parent,
        layout.inputs_dir,
        layout.artifacts_dir,
        layout.logs_dir,
        layout.reports_dir,
        layout.collections_dir,
        layout.drafts_dir,
        layout.light_models_dir,
        layout.full_models_dir,
        layout.light_adapter_data_dir,
        layout.full_adapter_data_dir,
        layout.archive_dir.parent,
        layout.matlab_scripts_dir,
        layout.matlab_helpers_dir,
        layout.matlab_adapters_dir,
    ):
        if dry_run:
            continue
        path.mkdir(parents=True, exist_ok=True)
    return layout


def log_message(logger: Logger, message: str) -> None:
    if logger is not None:
        logger(message)


def _format_command(command: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in command)


def _run_command(command: list[str], cwd: Path, logger: Logger, dry_run: bool = False) -> None:
    log_message(logger, f"$ {_format_command(command)}")
    if dry_run:
        return

    process = subprocess.Popen(
        command,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert process.stdout is not None
    for line in process.stdout:
        log_message(logger, line.rstrip("\n"))
    process.wait()
    if process.returncode != 0:
        raise RuntimeSetupError(
            f"Command failed with exit code {process.returncode}: {_format_command(command)}"
        )


def ensure_gecko_checkout(
    workspace: Path,
    gecko_dir: Path,
    repo_url: str = DEFAULT_GECKO_REPO,
    ref: str | None = None,
    logger: Logger = None,
    dry_run: bool = False,
) -> Path:
    ensure_runtime_layout(workspace, dry_run=dry_run)

    if gecko_dir.is_dir():
        if (gecko_dir / ".git").exists():
            log_message(logger, f"Using existing GECKO checkout: {gecko_dir}")
            return gecko_dir
        if (gecko_dir / "src").is_dir() and (gecko_dir / "README.md").is_file():
            log_message(logger, f"Using existing GECKO directory: {gecko_dir}")
            return gecko_dir
        raise RuntimeSetupError(
            f"GECKO directory exists but does not look like a usable checkout: {gecko_dir}"
        )

    git_exe = shutil.which("git")
    if not git_exe:
        raise RuntimeSetupError("git was not found on PATH, so GECKO cannot be cloned automatically.")

    clone_command = [git_exe, "clone", "--depth", "1"]
    if ref:
        clone_command.extend(["--branch", ref])
    clone_command.extend([repo_url, str(gecko_dir)])
    log_message(logger, f"Cloning GECKO into {gecko_dir}")
    _run_command(clone_command, cwd=workspace, logger=logger, dry_run=dry_run)
    return gecko_dir


def build_runtime_env(workspace: Path, gecko_dir: Path) -> dict[str, str]:
    layout = get_workspace_layout(workspace)
    env = os.environ.copy()
    env["MODELYML_WORKSPACE_ROOT"] = str(layout.workspace)
    env["MODELYML_INPUTS_DIR"] = str(layout.inputs_dir)
    env["MODELYML_ARTIFACTS_DIR"] = str(layout.artifacts_dir)
    env["MODELYML_DRAFTS_DIR"] = str(layout.drafts_dir)
    env["MODELYML_LIGHT_MODELS_DIR"] = str(layout.light_models_dir)
    env["MODELYML_FULL_MODELS_DIR"] = str(layout.full_models_dir)
    env["MODELYML_LIGHT_ADAPTER_DIR"] = str(layout.light_adapter_dir)
    env["MODELYML_FULL_ADAPTER_DIR"] = str(layout.full_adapter_dir)
    env["MODELYML_REPORT_DIR"] = str(layout.reports_dir)
    env["GECKO_MAIN_DIR"] = str(gecko_dir)
    return env


def legacy_generated_paths(workspace: Path) -> tuple[Path, ...]:
    layout = get_workspace_layout(workspace)
    return (
        workspace / "gecko_adapter",
        workspace / "gecko_full_adapter",
        workspace / "DATA_INPUT_OUTPUT",
        workspace / "pipeline_logs",
        workspace / "genome_for_carveme.faa",
        workspace / "genome_for_carveme_header_map.tsv",
        workspace / "model.xml",
        workspace / "model.yml",
        workspace / "model_nogapfill.xml",
        workspace / "rxn_to_ec.csv",
        workspace / "ecModel.yml",
        workspace / "ecModel_light.yml",
        workspace / "ecModel_kcat.yml",
        workspace / "ecModel_full.yml",
        workspace / "ecModel_full_kcat.yml",
        layout.workspace / "cmd",
    )
