import argparse
import csv
import hashlib
import json
from pathlib import Path
import ssl
import sys
from typing import Any
import urllib.error
import urllib.parse
import urllib.request

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from modelyml.pipeline_runtime import get_workspace_layout, infer_workspace_root  # type: ignore[no-redef]


UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
DEFAULT_CACHE_PATH = Path("metaMet/data/processed/predictkcat/uniprot_ec_cache.json")
DEFAULT_TIMEOUT_SECONDS = 20
DEFAULT_PAGE_SIZE = 25
DEFAULT_MAX_PAGES = 4
SSL_FALLBACK_USED = False


def normalize_header_value(value: str) -> str:
    return "_".join(str(value).strip().split())


def normalize_lookup_value(value: str) -> str:
    return " ".join(str(value).strip().lower().split())


def normalize_sequence(value: str) -> str:
    return str(value).strip().replace(" ", "").replace("\n", "").upper()


def auto_detect_column(columns: list[str], candidates: tuple[str, ...]) -> str | None:
    normalized = {c.lower().replace(" ", "").replace("_", ""): c for c in columns}
    for candidate in candidates:
        if candidate in normalized:
            return normalized[candidate]
    return None


def normalize_cli_args(argv: list[str]) -> list[str]:
    legacy_flags = {
        "input": "--input",
        "output": "--output",
        "organism": "--organism",
        "target-organism": "--target-organism",
        "target_organism": "--target-organism",
        "organism-column": "--organism-column",
        "organism_column": "--organism-column",
        "write-enriched-csv": "--write-enriched-csv",
        "write_enriched_csv": "--write-enriched-csv",
        "cache": "--cache",
        "skip-api-organism": "--skip-api-organism",
        "skip_api_organism": "--skip-api-organism",
    }

    normalized: list[str] = []
    for arg in argv:
        if arg.startswith("--"):
            normalized.append(arg)
            continue
        normalized.append(legacy_flags.get(arg, arg))
    return normalized


def load_cache(cache_path: Path) -> dict[str, dict[str, str]]:
    if not cache_path.exists():
        return {}
    try:
        data = json.loads(cache_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    if not isinstance(data, dict):
        return {}
    return {str(key): value for key, value in data.items() if isinstance(value, dict)}


def save_cache(cache_path: Path, cache: dict[str, dict[str, str]]) -> None:
    filtered_cache = {
        key: value
        for key, value in cache.items()
        if value.get("Organism") or value.get("Resolved Protein Sequence")
    }
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(filtered_cache, indent=2, sort_keys=True), encoding="utf-8")


def uniprot_search(query: str, *, size: int, offset: int = 0) -> list[dict[str, Any]]:
    global SSL_FALLBACK_USED

    params = urllib.parse.urlencode(
        {
            "query": query,
            "format": "json",
            "size": str(size),
            "offset": str(offset),
        }
    )
    request = urllib.request.Request(
        f"{UNIPROT_SEARCH_URL}?{params}",
        headers={"User-Agent": "ModelYML/Create_FASTA.py"},
    )
    try:
        with urllib.request.urlopen(request, timeout=DEFAULT_TIMEOUT_SECONDS) as response:
            payload = json.load(response)
    except urllib.error.URLError as error:
        reason = getattr(error, "reason", None)
        if isinstance(reason, ssl.SSLCertVerificationError):
            SSL_FALLBACK_USED = True
            insecure_context = ssl._create_unverified_context()
            with urllib.request.urlopen(request, timeout=DEFAULT_TIMEOUT_SECONDS, context=insecure_context) as response:
                payload = json.load(response)
        else:
            raise
    results = payload.get("results", [])
    return results if isinstance(results, list) else []


def build_uniprot_query(ec_number: str, organism: str | None = None, *, reviewed_only: bool) -> str:
    parts = [f"ec:{ec_number}"]
    if organism:
        parts.append(f'organism_name:"{organism}"')
    if reviewed_only:
        parts.append("reviewed:true")
    return " AND ".join(parts)


def empty_resolution(sequence_source: str = "unresolved") -> dict[str, str]:
    return {
        "Organism": "",
        "UniProt Accession": "",
        "Sequence Source": sequence_source,
        "Sequence Match": "unresolved",
        "Resolved Protein Sequence": "",
    }


def resolution_from_entry(entry: dict[str, Any], *, sequence_source: str, sequence_match: str) -> dict[str, str]:
    organism = entry.get("organism", {}) or {}
    sequence = entry.get("sequence", {}) or {}
    return {
        "Organism": str(organism.get("scientificName", "") or "").strip(),
        "UniProt Accession": str(entry.get("primaryAccession", "") or "").strip(),
        "Sequence Source": sequence_source,
        "Sequence Match": sequence_match,
        "Resolved Protein Sequence": normalize_sequence(sequence.get("value", "")),
    }


def resolve_metadata_for_existing_sequence(
    ec_number: str,
    sequence: str,
    cache: dict[str, dict[str, str]],
) -> dict[str, str]:
    cache_key = f"existing::{ec_number}::{hashlib.md5(sequence.encode('utf-8')).hexdigest()}"
    cached = cache.get(cache_key)
    if cached and (cached.get("Organism") or cached.get("Resolved Protein Sequence")):
        return cached

    if not ec_number:
        return empty_resolution()

    for reviewed_only in (True, False):
        query = build_uniprot_query(ec_number, reviewed_only=reviewed_only)
        for page in range(DEFAULT_MAX_PAGES):
            try:
                entries = uniprot_search(query, size=DEFAULT_PAGE_SIZE, offset=page * DEFAULT_PAGE_SIZE)
            except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
                entries = []
            if not entries:
                break
            for entry in entries:
                if normalize_sequence((entry.get("sequence", {}) or {}).get("value", "")) == sequence:
                    resolved = resolution_from_entry(
                        entry,
                        sequence_source="uniprot_ec_exact_sequence",
                        sequence_match="exact",
                    )
                    cache[cache_key] = resolved
                    return resolved

    for reviewed_only in (True, False):
        query = build_uniprot_query(ec_number, reviewed_only=reviewed_only)
        try:
            entries = uniprot_search(query, size=1)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            entries = []
        if entries:
            resolved = resolution_from_entry(
                entries[0],
                sequence_source="uniprot_ec_best_effort",
                sequence_match="best_effort",
            )
            cache[cache_key] = resolved
            return resolved

    return empty_resolution()


def resolve_target_organism_sequence(
    ec_number: str,
    organism: str,
    cache: dict[str, dict[str, str]],
) -> dict[str, str]:
    cache_key = f"target::{ec_number}::{normalize_lookup_value(organism)}"
    cached = cache.get(cache_key)
    if cached and cached.get("Resolved Protein Sequence"):
        return cached

    if not ec_number:
        return empty_resolution()

    for reviewed_only in (True, False):
        query = build_uniprot_query(ec_number, organism, reviewed_only=reviewed_only)
        try:
            entries = uniprot_search(query, size=1)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            entries = []
        if entries:
            resolved = resolution_from_entry(
                entries[0],
                sequence_source="uniprot_target_organism_reviewed" if reviewed_only else "uniprot_target_organism",
                sequence_match="target_organism",
            )
            cache[cache_key] = resolved
            return resolved

    return empty_resolution()


def main() -> None:
    layout = get_workspace_layout(infer_workspace_root(Path(__file__)))
    parser = argparse.ArgumentParser(
        description=(
            "Create genome.faa from the predictkcat CSV, automatically carry real organism metadata into FASTA headers, "
            "and optionally fetch organism-specific sequences from UniProt for a single-organism model."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("metaMet/data/processed/predictkcat/input_with_ec.csv"),
        help="Input CSV file.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=layout.genome_faa,
        help="Output FASTA file.",
    )
    parser.add_argument(
        "--organism",
        help="Single organism name to add to every FASTA header, e.g. 'Escherichia coli'.",
    )
    parser.add_argument(
        "--target-organism",
        help=(
            "Fetch organism-specific protein sequences from UniProt using the EC number and build a real single-organism FASTA. "
            "This is the recommended mode when you want model.xml/model.yml for one organism."
        ),
    )
    parser.add_argument(
        "--organism-column",
        help="CSV column containing organism names. If omitted, a common organism column is auto-detected.",
    )
    parser.add_argument(
        "--write-enriched-csv",
        type=Path,
        help="Optional CSV path to write an enriched table with Organism, UniProt Accession, and Sequence Source columns.",
    )
    parser.add_argument(
        "--cache",
        type=Path,
        default=DEFAULT_CACHE_PATH,
        help="JSON cache for UniProt EC lookups.",
    )
    parser.add_argument(
        "--skip-api-organism",
        action="store_true",
        help="Skip UniProt organism lookup when the input CSV does not already contain organism information.",
    )
    args = parser.parse_args(normalize_cli_args(sys.argv[1:]))

    with args.input.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames or []
        rows = list(reader)

    seq_col_candidates = [c for c in columns if "protein" in c.lower() and "sequence" in c.lower()]
    if not seq_col_candidates:
        raise ValueError(f"Couldn't find protein sequence column. Columns: {columns}")
    seq_col = seq_col_candidates[0]

    ec_col = auto_detect_column(columns, ("ecnumber", "ec"))

    organism_col = args.organism_column
    if organism_col is None:
        organism_col = auto_detect_column(
            columns,
            (
                "organism",
                "organismname",
                "species",
                "scientificname",
                "taxonname",
                "strain",
            ),
        )
    if organism_col and organism_col not in columns:
        raise ValueError(f"Organism column '{organism_col}' not found. Columns: {columns}")

    cache = load_cache(args.cache)
    seen = set()
    records = []
    enriched_rows: list[dict[str, str]] = []
    api_queries = 0
    organism_hits = 0
    target_hits = 0
    skipped_for_target = 0

    for i, row in enumerate(rows):
        seq = normalize_sequence(row.get(seq_col, ""))
        if not seq or seq == "nan":
            continue

        ec = str(row.get(ec_col, "")).strip() if ec_col else "NA"
        if not ec or ec.lower() == "nan":
            ec = "NA"

        resolved = empty_resolution(sequence_source="input_csv")
        if args.target_organism:
            api_queries += 1
            resolved = resolve_target_organism_sequence(ec, args.target_organism, cache)
            resolved_sequence = normalize_sequence(resolved.get("Resolved Protein Sequence", ""))
            if not resolved_sequence:
                skipped_for_target += 1
                continue
            seq = resolved_sequence
            target_hits += 1
        else:
            organism_in_row = str(row.get(organism_col, "")).strip() if organism_col else ""
            if organism_in_row and organism_in_row.lower() != "nan":
                resolved["Organism"] = organism_in_row
                resolved["Sequence Source"] = "input_csv"
                resolved["Sequence Match"] = "input_column"
            elif not args.skip_api_organism and ec != "NA":
                api_queries += 1
                resolved = resolve_metadata_for_existing_sequence(ec, seq, cache)
                if resolved.get("Organism"):
                    organism_hits += 1

        organism_value = args.organism
        resolved_organism = str(resolved.get("Organism", "")).strip()
        if resolved_organism:
            organism_value = resolved_organism
        elif args.target_organism:
            organism_value = args.target_organism

        enriched_row = {key: str(value) for key, value in row.items()}
        enriched_row[seq_col] = seq
        enriched_row["Organism"] = organism_value or ""
        enriched_row["UniProt Accession"] = resolved.get("UniProt Accession", "")
        enriched_row["Sequence Source"] = resolved.get("Sequence Source", "")
        enriched_row["Sequence Match"] = resolved.get("Sequence Match", "")
        enriched_rows.append(enriched_row)

        if seq in seen:
            continue
        seen.add(seq)

        header_parts = [f">seq{i+1}"]
        if organism_value:
            header_parts.append(f"ORG={normalize_header_value(organism_value)}")

        header_parts.append(f"EC={ec}")
        accession = str(resolved.get("UniProt Accession", "")).strip()
        if accession:
            header_parts.append(f"UP={accession}")

        header = "|".join(header_parts)
        records.append(header + "\n" + seq + "\n")

    if args.write_enriched_csv:
        output_columns = list(columns)
        for extra_column in ("Organism", "UniProt Accession", "Sequence Source", "Sequence Match"):
            if extra_column not in output_columns:
                output_columns.append(extra_column)
        args.write_enriched_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.write_enriched_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=output_columns)
            writer.writeheader()
            writer.writerows(enriched_rows)

    save_cache(args.cache, cache)
    args.output.write_text("".join(records))
    if args.write_enriched_csv:
        print(f"Wrote enriched CSV to {args.write_enriched_csv}")
    if SSL_FALLBACK_USED:
        print("Warning: UniProt lookups used an SSL fallback because the local Python certificate store is incomplete.")
    if args.target_organism:
        print(
            f"Resolved {target_hits} rows for target organism '{args.target_organism}' and skipped {skipped_for_target} rows with no UniProt match."
        )
    elif api_queries:
        print(f"Resolved organism metadata from UniProt for {organism_hits} rows.")
    print(f"Wrote {len(records)} unique protein sequences to {args.output}")


if __name__ == "__main__":
    main()