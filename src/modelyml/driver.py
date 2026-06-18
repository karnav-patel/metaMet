from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from modelyml.config import PipelineConfig, load_pipeline_config
from modelyml.pipeline_runtime import get_workspace_layout, infer_workspace_root, legacy_generated_paths


class DriverError(RuntimeError):
    pass


def parse_args() -> argparse.Namespace:
    workspace = infer_workspace_root(Path(__file__), levels_up=3)
    parser = argparse.ArgumentParser(description='Run the organized ModelYML workflows from a single config file.')
    parser.add_argument('--config', type=Path, default=workspace / 'config' / 'pipeline.toml', help='TOML config file.')
    return parser.parse_args()


def _append_optional(command: list[str], flag: str, value: str | None) -> None:
    if value:
        command.extend([flag, value])


def _run(command: list[str], workspace: Path) -> None:
    completed = subprocess.run(command, cwd=workspace)
    if completed.returncode != 0:
        raise DriverError(f"Command failed with exit code {completed.returncode}: {' '.join(command)}")


def cleanup_legacy_root_outputs(workspace: Path, dry_run: bool) -> None:
    for path in legacy_generated_paths(workspace):
        if not path.exists():
            continue
        print(f"Remove legacy root artifact: {path}")
        if dry_run:
            continue
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        else:
            path.unlink()


def build_base_command(script_path: Path, config: PipelineConfig, workspace: Path) -> list[str]:
    python_cmd = config.tools.python_cmd or sys.executable
    command = [python_cmd, str(script_path), '--workspace', str(workspace)]
    command.extend(['--gecko-dir', str(config.gecko.directory), '--gecko-repo', config.gecko.repo])
    if config.gecko.ref:
        command.extend(['--gecko-ref', config.gecko.ref])
    _append_optional(command, '--matlab-cmd', config.tools.matlab_cmd)
    _append_optional(command, '--carve-cmd', config.tools.carve_cmd)
    _append_optional(command, '--matlab-mode', config.tools.matlab_mode)
    if config.run.resume:
        command.append('--resume')
    if config.run.dry_run:
        command.append('--dry-run')
    return command


def run_all_modes(config: PipelineConfig, workspace: Path) -> None:
    layout = get_workspace_layout(workspace)
    script_path = layout.python_pipelines_dir / 'collect_outputs.py'
    command = build_base_command(script_path, config, workspace)
    command.extend(['--output-root', str((workspace / config.paths.collections_dir).resolve())])
    command.extend(['--gapfill-media', config.run.gapfill_media])
    _run(command, workspace)


def run_single_mode(config: PipelineConfig, workspace: Path, mode: str) -> None:
    layout = get_workspace_layout(workspace)
    if mode == 'light':
        script_path = layout.light_pipeline_script
    elif mode == 'full':
        script_path = layout.full_pipeline_script
    else:
        raise DriverError(f'Unsupported run mode: {mode}')

    base = build_base_command(script_path, config, workspace)
    variants: list[list[str]] = []
    if config.run.include_no_gapfill:
        variants.append(base.copy())
    if config.run.include_gapfill:
        gapfill_cmd = base.copy()
        gapfill_cmd.extend(['--gapfill-media', config.run.gapfill_media])
        variants.append(gapfill_cmd)

    if not variants:
        raise DriverError('Nothing to run: both include_gapfill and include_no_gapfill are disabled.')

    for command in variants:
        _run(command, workspace)


def main() -> int:
    args = parse_args()
    config_path = args.config.resolve()
    workspace = infer_workspace_root(Path(__file__), levels_up=3)
    config = load_pipeline_config(config_path)

    if config.run.cleanup_legacy_root:
        cleanup_legacy_root_outputs(workspace, config.run.dry_run)

    mode = config.run.mode
    if mode == 'all':
        run_all_modes(config, workspace)
    else:
        run_single_mode(config, workspace, mode)
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except DriverError as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        raise SystemExit(1)
