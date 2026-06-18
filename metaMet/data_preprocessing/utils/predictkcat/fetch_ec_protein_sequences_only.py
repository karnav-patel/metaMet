# fetch_ec_protein_sequences_only.py
import sys
import csv
import time
import requests
from pathlib import Path

def get_project_root():
    current_file_path = Path(__file__).resolve()
    for parent in current_file_path.parents:
        if parent.name == "metaMet":
            return parent
    raise RuntimeError("Project root folder 'metaMet' not found.")

project_root = get_project_root()
data_preproc_dir = project_root / "data_preprocessing"
sys.path.insert(0, str(data_preproc_dir))

import config.config as config

def get_ecs_missing_kcat():
    brenda_file = config.brenda_kcat_csv
    metacyc_file = config.metacyc_kcat_csv
    sabio_file  = config.sabio_kcat_csv
    merged_file = config.merged_extracted_ids_output_path

    all_ecs = set()
    with merged_file.open("r", encoding="utf-8") as f:
        for line in f:
            ec = line.strip()
            if ec:
                all_ecs.add(ec)

    ecs_with_kcat = set()

    def iter_rows(path: Path, delimiter=","):
        with path.open("r", encoding="utf-8", newline="") as f:
            yield from csv.DictReader(f, delimiter=delimiter)

    for row in iter_rows(brenda_file, delimiter=","):
        ec   = row.get("ec_number", "").strip()
        kcat = row.get("kcat", "").strip()
        if ec and kcat:
            ecs_with_kcat.add(ec)

    for row in iter_rows(metacyc_file, delimiter=","):
        ec   = row.get("ec_number", "").strip()
        kcat = row.get("kcat", "").strip()
        if ec and kcat:
            ecs_with_kcat.add(ec)

    for row in iter_rows(sabio_file, delimiter="\t"):
        if row.get("Type", "").strip().lower() == "kcat":
            ec = row.get("ECNumber", "").strip()
            if ec:
                ecs_with_kcat.add(ec)

    missing_ecs = sorted(all_ecs - ecs_with_kcat)
    return missing_ecs

def uniprot_search_ec(ec):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=ec:{ec}&fields=accession&format=json&size=1"
    j = requests.get(url, timeout=5).json().get("results",[])
    return j[0]["primaryAccession"] if j else None

def fetch_uniprot_sequence(acc):
    url1 = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    try:
        return requests.get(url1, timeout=5).json()["sequence"]["value"]
    except: pass
    url2 = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    try:
        text = requests.get(url2, timeout=5).text.splitlines()[1:]
        return "".join(text)
    except: pass
    return ""

def get_kegg_seq(ec):
    link = f"http://rest.kegg.jp/link/uniprot/enzyme:{ec}"
    try:
        txt  = requests.get(link, timeout=5).text.strip().splitlines()
    except:
        return ""
    if not txt:
        return ""
    up = txt[0].split("\t")[1].replace("up:","")
    get = f"http://rest.kegg.jp/get/up:{up}"
    seq, on = "", False
    try:
        for L in requests.get(get, timeout=5).text.splitlines():
            if L.startswith("AASEQ"):
                on=True; seq+=L.split()[1]
            elif on:
                if L.startswith(" "): seq+=L.strip()
                else: break
    except:
        return ""
    return seq

def fetch_protein_sequence(ec):
    acc = uniprot_search_ec(ec)
    if acc:
        seq = fetch_uniprot_sequence(acc)
        if seq: return seq
    return get_kegg_seq(ec)

def main():
    ecs = [e.strip().lower() for e in get_ecs_missing_kcat()]
    out_csv = config.predictkcat_with_ec_protein_sequences  # CSV with columns: ec_number,protein_sequence

    rows = []
    for ec in ecs:
        seq = fetch_protein_sequence(ec)
        rows.append({"ec_number": ec, "protein_sequence": seq})
        time.sleep(0.2)

    # Write CSV
    with open(out_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["ec_number","protein_sequence"])
        w.writeheader()
        for r in rows:
            w.writerow(r)

if __name__ == "__main__":
    main()
