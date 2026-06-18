from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from modelyml.pipeline_runtime import get_workspace_layout, infer_workspace_root  # type: ignore[no-redef]

_DEFAULT_LAYOUT = get_workspace_layout(infer_workspace_root(Path(__file__)))
DEFAULT_INPUT = _DEFAULT_LAYOUT.genome_faa
DEFAULT_OUTPUT = _DEFAULT_LAYOUT.sanitized_fasta
DEFAULT_MAP = _DEFAULT_LAYOUT.sanitized_fasta_map


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rewrite FASTA headers into CarveMe-safe identifiers. "
            "This avoids parser failures caused by headers containing pipes, parentheses, spaces, or boolean-like text."
        )
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Input FASTA path.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Sanitized FASTA output path.")
    parser.add_argument(
        "--map-out",
        type=Path,
        default=DEFAULT_MAP,
        help="TSV file recording original headers and sanitized IDs.",
    )
    return parser.parse_args()


def iter_fasta_records(text: str):
    header: str | None = None
    seq_lines: list[str] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, seq_lines
            header = line[1:]
            seq_lines = []
        else:
            if header is None:
                continue
            seq_lines.append(line)

    if header is not None:
        yield header, seq_lines


def make_safe_identifier(header: str, index: int, used: set[str]) -> str:
    first_token = header.split()[0]
    base = first_token.split("|")[0]
    base = re.sub(r"[^A-Za-z0-9_]", "_", base)
    base = re.sub(r"_+", "_", base).strip("_")

    if not base:
        base = f"seq_{index}"
    if not re.match(r"[A-Za-z_]", base):
        base = f"seq_{base}"

    candidate = base
    suffix = 2
    while candidate in used:
        candidate = f"{base}_{suffix}"
        suffix += 1
    used.add(candidate)
    return candidate


def main() -> int:
    args = parse_args()
    input_path = args.input.resolve()
    output_path = args.output.resolve()
    map_path = args.map_out.resolve()

    if not input_path.is_file():
        print(f"ERROR: FASTA file not found: {input_path}", file=sys.stderr)
        return 1

    records = list(iter_fasta_records(input_path.read_text(encoding="utf-8")))
    if not records:
        print(f"ERROR: No FASTA records found in {input_path}", file=sys.stderr)
        return 1

    used: set[str] = set()
    output_lines: list[str] = []
    rows: list[tuple[str, str]] = []

    for index, (header, seq_lines) in enumerate(records, start=1):
        safe_id = make_safe_identifier(header, index, used)
        output_lines.append(f">{safe_id}")
        output_lines.extend(seq_lines)
        rows.append((safe_id, header))

    output_path.write_text("\n".join(output_lines) + "\n", encoding="utf-8")
    map_path.parent.mkdir(parents=True, exist_ok=True)
    with map_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["safe_id", "original_header"])
        writer.writerows(rows)

    print(f"Input file: {input_path}")
    print(f"Output file: {output_path}")
    print(f"Header map: {map_path}")
    print(f"Records rewritten: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
