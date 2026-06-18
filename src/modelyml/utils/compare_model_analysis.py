from __future__ import annotations

import argparse
import difflib
import math
import re
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from modelyml.pipeline_runtime import (  # type: ignore[no-redef]
    DEFAULT_ARCHIVE_DIR,
    get_workspace_layout,
    infer_workspace_root,
)


TOP_LEVEL_SECTION_RE = re.compile(r"^-\s+([A-Za-z0-9_]+):\s*$")
ID_LINE_RE = re.compile(r'^\s+-\s+id:\s+"?([^"\n]+)"?\s*$')
KCAT_LINE_RE = re.compile(r"^\s+-\s+kcat:\s+([^\s#]+)")


@dataclass
class SectionDiff:
    name: str
    current_count: int
    old_count: int
    added: list[str]
    removed: list[str]


@dataclass
class KcatStats:
    total: int
    nonzero: int
    zeros: int
    minimum: float | None
    maximum: float | None
    mean: float | None


def parse_args() -> argparse.Namespace:
    root = infer_workspace_root(Path(__file__))
    default_report = get_workspace_layout(root).reports_dir / "compare_analysis_report.md"
    parser = argparse.ArgumentParser(
        description="Compare model files and write a Markdown report."
    )
    parser.add_argument("--workspace", type=Path, default=root, help="Workspace root.")
    parser.add_argument(
        "--report",
        type=Path,
        default=default_report,
        help="Markdown report output path.",
    )
    parser.add_argument(
        "--max-items",
        type=int,
        default=25,
        help="Maximum added/removed IDs to show per section.",
    )
    parser.add_argument(
        "--current-file",
        type=Path,
        default=None,
        help="Optional current file to compare instead of using the default archive/oldfiles report mode.",
    )
    parser.add_argument(
        "--reference-file",
        type=Path,
        default=None,
        help="Optional reference file paired with --current-file.",
    )
    parser.add_argument(
        "--kind",
        choices=("auto", "yaml", "xml"),
        default="auto",
        help="Comparison kind for single-pair mode. Defaults to auto-detection from the file suffix.",
    )
    parser.add_argument(
        "--report-title",
        default=None,
        help="Optional Markdown title for the generated report.",
    )
    parser.add_argument(
        "--report-description",
        default=None,
        help="Optional free-text description shown near the top of the report.",
    )
    return parser.parse_args()


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def parse_yaml_sections(path: Path) -> dict[str, list[str]]:
    sections: dict[str, list[str]] = {}
    current_section: str | None = None

    for line in read_text(path).splitlines():
        top_match = TOP_LEVEL_SECTION_RE.match(line)
        if top_match:
            current_section = top_match.group(1)
            sections.setdefault(current_section, [])
            continue

        if current_section is None:
            continue

        id_match = ID_LINE_RE.match(line)
        if id_match:
            sections[current_section].append(id_match.group(1))

    return sections


def parse_kcat_stats(path: Path) -> KcatStats:
    values: list[float] = []
    for line in read_text(path).splitlines():
        match = KCAT_LINE_RE.match(line)
        if not match:
            continue
        raw = match.group(1).strip().strip('"')
        try:
            value = float(raw)
        except ValueError:
            continue
        values.append(value)

    if not values:
        return KcatStats(0, 0, 0, None, None, None)

    nonzero_values = [value for value in values if value != 0]
    return KcatStats(
        total=len(values),
        nonzero=len(nonzero_values),
        zeros=len(values) - len(nonzero_values),
        minimum=min(values),
        maximum=max(values),
        mean=sum(values) / len(values),
    )


def compare_id_lists(name: str, current_ids: Iterable[str], old_ids: Iterable[str]) -> SectionDiff:
    current_set = set(current_ids)
    old_set = set(old_ids)
    return SectionDiff(
        name=name,
        current_count=len(current_set),
        old_count=len(old_set),
        added=sorted(current_set - old_set),
        removed=sorted(old_set - current_set),
    )


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def parse_xml_summary(path: Path) -> dict[str, list[str]]:
    root = ET.parse(path).getroot()
    summary = {"species": [], "reactions": [], "geneProducts": []}

    for element in root.iter():
        name = local_name(element.tag)
        if name == "species":
            identifier = element.attrib.get("id")
            if identifier:
                summary["species"].append(identifier)
        elif name == "reaction":
            identifier = element.attrib.get("id")
            if identifier:
                summary["reactions"].append(identifier)
        elif name == "geneProduct":
            identifier = element.attrib.get("id") or element.attrib.get("label")
            if identifier:
                summary["geneProducts"].append(identifier)

    return summary


def summarize_text_diff(current_path: Path, old_path: Path, max_lines: int = 80) -> list[str]:
    diff = list(
        difflib.unified_diff(
            read_text(old_path).splitlines(),
            read_text(current_path).splitlines(),
            fromfile=str(old_path.name),
            tofile=str(current_path.name),
            n=1,
            lineterm="",
        )
    )
    return diff[:max_lines]


def format_float(value: float | None) -> str:
    if value is None:
        return "n/a"
    if math.isclose(value, round(value), rel_tol=0, abs_tol=1e-12):
        return str(int(round(value)))
    return f"{value:.6g}"


def markdown_list(items: list[str], max_items: int) -> str:
    if not items:
        return "none"
    shown = items[:max_items]
    text = ", ".join(shown)
    if len(items) > max_items:
        text += f", ... (+{len(items) - max_items} more)"
    return text


def compare_yaml_pair(current_path: Path, reference_path: Path, max_items: int) -> list[str]:
    current_sections = parse_yaml_sections(current_path)
    old_sections = parse_yaml_sections(reference_path)

    lines = [f"## {current_path.name}", ""]
    lines.append(f"Current file: {current_path}")
    lines.append(f"Reference file: {reference_path}")
    lines.append("")
    lines.append("| Section | Current | Reference | Δ |")
    lines.append("|---|---:|---:|---:|")

    tracked_sections = [
        section
        for section in ["metabolites", "reactions", "genes"]
        if section in current_sections or section in old_sections
    ]

    diffs: list[SectionDiff] = []
    for section in tracked_sections:
        diff = compare_id_lists(section, current_sections.get(section, []), old_sections.get(section, []))
        diffs.append(diff)
        delta = diff.current_count - diff.old_count
        lines.append(f"| {section} | {diff.current_count} | {diff.old_count} | {delta:+d} |")

    if current_path.name == "ecModel_kcat.yml":
        current_kcat = parse_kcat_stats(current_path)
        old_kcat = parse_kcat_stats(reference_path)
        lines.append("")
        lines.append("### kcat summary")
        lines.append("")
        lines.append("| Metric | Current | Reference |")
        lines.append("|---|---:|---:|")
        lines.append(f"| Total kcat entries | {current_kcat.total} | {old_kcat.total} |")
        lines.append(f"| Nonzero kcats | {current_kcat.nonzero} | {old_kcat.nonzero} |")
        lines.append(f"| Zero kcats | {current_kcat.zeros} | {old_kcat.zeros} |")
        lines.append(f"| Min kcat | {format_float(current_kcat.minimum)} | {format_float(old_kcat.minimum)} |")
        lines.append(f"| Max kcat | {format_float(current_kcat.maximum)} | {format_float(old_kcat.maximum)} |")
        lines.append(f"| Mean kcat | {format_float(current_kcat.mean)} | {format_float(old_kcat.mean)} |")

    lines.append("")
    lines.append("### Added/removed IDs")
    lines.append("")
    for diff in diffs:
        lines.append(f"- **{diff.name}** added: {markdown_list(diff.added, max_items)}")
        lines.append(f"- **{diff.name}** removed: {markdown_list(diff.removed, max_items)}")
    lines.append("")
    lines.append("### Text diff excerpt")
    lines.append("")
    lines.append("```")
    lines.extend(summarize_text_diff(current_path, reference_path))
    lines.append("```")
    lines.append("")
    return lines


def compare_xml_pair(current_path: Path, reference_path: Path, max_items: int) -> list[str]:
    current_summary = parse_xml_summary(current_path)
    old_summary = parse_xml_summary(reference_path)

    lines = [f"## {current_path.name}", ""]
    lines.append(f"Current file: {current_path}")
    lines.append(f"Reference file: {reference_path}")
    lines.append("")
    lines.append("| Section | Current | Reference | Δ |")
    lines.append("|---|---:|---:|---:|")

    diffs: list[SectionDiff] = []
    for section in ["species", "reactions", "geneProducts"]:
        diff = compare_id_lists(section, current_summary.get(section, []), old_summary.get(section, []))
        diffs.append(diff)
        delta = diff.current_count - diff.old_count
        lines.append(f"| {section} | {diff.current_count} | {diff.old_count} | {delta:+d} |")

    lines.append("")
    lines.append("### Added/removed IDs")
    lines.append("")
    for diff in diffs:
        lines.append(f"- **{diff.name}** added: {markdown_list(diff.added, max_items)}")
        lines.append(f"- **{diff.name}** removed: {markdown_list(diff.removed, max_items)}")
    lines.append("")
    lines.append("### Text diff excerpt")
    lines.append("")
    lines.append("```")
    lines.extend(summarize_text_diff(current_path, reference_path))
    lines.append("```")
    lines.append("")
    return lines


def resolve_path(path: Path, workspace: Path) -> Path:
    return path if path.is_absolute() else workspace / path


def detect_kind(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".xml":
        return "xml"
    if suffix in {".yml", ".yaml"}:
        return "yaml"
    raise ValueError(f"Could not infer comparison kind from file suffix: {path}")


def main() -> int:
    args = parse_args()
    workspace = args.workspace.resolve()
    report_path = resolve_path(args.report, workspace).resolve()
    layout = get_workspace_layout(workspace)
    old_dir = workspace / DEFAULT_ARCHIVE_DIR

    if bool(args.current_file) != bool(args.reference_file):
        raise SystemExit("ERROR: --current-file and --reference-file must be provided together.")

    if args.current_file and args.reference_file:
        current_path = resolve_path(args.current_file, workspace).resolve()
        reference_path = resolve_path(args.reference_file, workspace).resolve()
        kind = args.kind if args.kind != "auto" else detect_kind(current_path)
        report_path.parent.mkdir(parents=True, exist_ok=True)

        lines = [
            f"# {args.report_title or 'Comparison report'}",
            "",
            f"Workspace: {workspace}",
            "",
        ]
        if args.report_description:
            lines.append(args.report_description)
            lines.append("")

        if not current_path.exists() or not reference_path.exists():
            missing = [str(path) for path in (current_path, reference_path) if not path.exists()]
            raise SystemExit("ERROR: Missing comparison file(s): " + ", ".join(missing))

        if kind == "yaml":
            lines.extend(compare_yaml_pair(current_path, reference_path, args.max_items))
        else:
            lines.extend(compare_xml_pair(current_path, reference_path, args.max_items))

        report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"Wrote report: {report_path}")
        return 0

    comparisons = [
        (layout.light_ec_model_kcat, old_dir / "ecModel_kcat.yml", "yaml"),
        (layout.draft_model_yaml, old_dir / "model.yml", "yaml"),
        (layout.draft_model_xml, old_dir / "model.xml", "xml"),
    ]

    optional = (layout.light_ec_model, old_dir / "ecModel.yml", "yaml")
    if optional[0].exists() and optional[1].exists():
        comparisons.append(optional)

    lines = [
        "# Compare analysis report",
        "",
        f"Workspace: {workspace}",
        f"Old reference folder: {old_dir}",
        "",
        f"This report compares the current one-organism run against the older multi-organism outputs stored in {old_dir.relative_to(workspace)}/.",
        "",
    ]

    for current_path, old_path, kind in comparisons:
        if not current_path.exists() or not old_path.exists():
            lines.append(f"## {current_path.name}")
            lines.append("")
            lines.append("Missing comparison file; skipped.")
            lines.append("")
            continue

        if kind == "yaml":
            lines.extend(compare_yaml_pair(current_path, old_path, args.max_items))
        else:
            lines.extend(compare_xml_pair(current_path, old_path, args.max_items))

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
