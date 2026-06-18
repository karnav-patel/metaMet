from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from modelyml.pipeline_runtime import get_workspace_layout, infer_workspace_root  # type: ignore[no-redef]

DEFAULT_ORGANISM = "ORG=Escherichia_coli_(strain_K12)"


def parse_args() -> argparse.Namespace:
    layout = get_workspace_layout(infer_workspace_root(Path(__file__)))
    parser = argparse.ArgumentParser(
        description=(
            "Filter a FASTA file so only records whose header contains the target organism remain. "
            "By default this edits genome.faa in place and keeps a .bak backup."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=layout.genome_faa,
        help="Input FASTA file. Defaults to genome.faa.",
    )
    parser.add_argument(
        "--organism",
        default=DEFAULT_ORGANISM,
        help="Exact header token to keep. Defaults to ORG=Escherichia_coli_(strain_K12).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output FASTA path. If omitted, the input file is overwritten.",
    )
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Do not create a .bak backup when overwriting the input file.",
    )
    return parser.parse_args()


def iter_fasta_records(text: str):
    header: str | None = None
    seq_lines: list[str] = []

    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")
        if line.startswith(">"):
            if header is not None:
                yield header, seq_lines
            header = line
            seq_lines = []
        else:
            if header is None:
                continue
            seq_lines.append(line)

    if header is not None:
        yield header, seq_lines


def main() -> int:
    args = parse_args()
    input_path = args.input.resolve()
    output_path = args.output.resolve() if args.output else input_path

    if not input_path.is_file():
        print(f"ERROR: Input FASTA not found: {input_path}", file=sys.stderr)
        return 1

    text = input_path.read_text(encoding="utf-8")
    records = list(iter_fasta_records(text))
    if not records:
        print(f"ERROR: No FASTA records found in {input_path}", file=sys.stderr)
        return 1

    kept: list[str] = []
    kept_count = 0
    removed_count = 0

    for header, seq_lines in records:
        if args.organism in header:
            kept.append(header)
            kept.extend(seq_lines)
            kept_count += 1
        else:
            removed_count += 1

    if kept_count == 0:
        print(
            f"ERROR: No records matched organism token '{args.organism}'. File was not changed.",
            file=sys.stderr,
        )
        return 1

    output_text = "\n".join(kept) + "\n"

    if output_path == input_path and not args.no_backup:
        backup_path = input_path.with_suffix(input_path.suffix + ".bak")
        shutil.copy2(input_path, backup_path)
        print(f"Backup created: {backup_path}")

    output_path.write_text(output_text, encoding="utf-8")

    print(f"Input file: {input_path}")
    print(f"Output file: {output_path}")
    print(f"Target organism: {args.organism}")
    print(f"Records kept: {kept_count}")
    print(f"Records removed: {removed_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
