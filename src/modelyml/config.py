from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tomllib

from modelyml.pipeline_runtime import DEFAULT_GECKO_DIR, DEFAULT_GECKO_REPO


@dataclass(frozen=True)
class RunOptions:
    mode: str
    include_gapfill: bool
    include_no_gapfill: bool
    resume: bool
    dry_run: bool
    cleanup_legacy_root: bool
    gapfill_media: str


@dataclass(frozen=True)
class ToolOptions:
    python_cmd: str | None
    matlab_cmd: str | None
    carve_cmd: str | None
    matlab_mode: str | None


@dataclass(frozen=True)
class GeckoOptions:
    directory: Path
    repo: str
    ref: str | None


@dataclass(frozen=True)
class PathOptions:
    collections_dir: Path


@dataclass(frozen=True)
class PipelineConfig:
    run: RunOptions
    tools: ToolOptions
    gecko: GeckoOptions
    paths: PathOptions


def _get_bool(section: dict, key: str, default: bool) -> bool:
    value = section.get(key, default)
    return bool(value)


def _get_optional_str(section: dict, key: str) -> str | None:
    value = section.get(key)
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def load_pipeline_config(config_path: Path) -> PipelineConfig:
    data = tomllib.loads(config_path.read_text(encoding='utf-8'))

    run_section = data.get('run', {})
    tools_section = data.get('tools', {})
    gecko_section = data.get('gecko', {})
    paths_section = data.get('paths', {})

    run = RunOptions(
        mode=str(run_section.get('mode', 'all')).strip().lower(),
        include_gapfill=_get_bool(run_section, 'include_gapfill', True),
        include_no_gapfill=_get_bool(run_section, 'include_no_gapfill', True),
        resume=_get_bool(run_section, 'resume', True),
        dry_run=_get_bool(run_section, 'dry_run', False),
        cleanup_legacy_root=_get_bool(run_section, 'cleanup_legacy_root', True),
        gapfill_media=str(run_section.get('gapfill_media', 'M9')).strip() or 'M9',
    )

    tools = ToolOptions(
        python_cmd=_get_optional_str(tools_section, 'python_cmd'),
        matlab_cmd=_get_optional_str(tools_section, 'matlab_cmd'),
        carve_cmd=_get_optional_str(tools_section, 'carve_cmd'),
        matlab_mode=_get_optional_str(tools_section, 'matlab_mode'),
    )

    gecko = GeckoOptions(
        directory=Path(str(gecko_section.get('dir', DEFAULT_GECKO_DIR))),
        repo=str(gecko_section.get('repo', DEFAULT_GECKO_REPO)).strip() or DEFAULT_GECKO_REPO,
        ref=_get_optional_str(gecko_section, 'ref'),
    )

    paths = PathOptions(
        collections_dir=Path(str(paths_section.get('collections_dir', 'artifacts/collections'))),
    )

    return PipelineConfig(run=run, tools=tools, gecko=gecko, paths=paths)
