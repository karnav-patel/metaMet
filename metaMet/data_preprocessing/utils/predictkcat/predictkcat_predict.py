#!/usr/bin/env python3
import os
import sys
import subprocess
import shutil
import zipfile
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

def clone_dlkcat():
    repo_url   = "https://github.com/SysBioChalmers/DLKcat.git"
    base_dir   = Path(config.predictkcat_DLKcat)       # …/predictkcat/DLKcat
    repo_dest  = base_dir / "DLKcat"                   # …/predictkcat/DLKcat/DLKcat

    # if the inner DLKcat already exists, we assume it's been cloned
    if repo_dest.exists():
        print(f"→ '{repo_dest}' already exists — skipping clone.")
        return

    # make sure the outer DLKcat folder exists
    if not base_dir.exists():
        print(f"→ Creating directory {base_dir}")
        base_dir.mkdir(parents=True, exist_ok=True)

    print(f"→ Cloning {repo_url} into {repo_dest}")
    proc = subprocess.Popen(
        ["git", "clone", repo_url, str(repo_dest)],
        stdout=sys.stdout,
        stderr=sys.stderr
    )
    ret = proc.wait()
    if ret != 0:
        print(f"→ Git clone failed with exit code {ret}")
        sys.exit(ret)

    print("→ Clone completed successfully.")


def replace_example_input():
    src = Path(config.predictkcat_input_tsv)
    if not src.exists():
        print(f"→ Source file {src} does not exist; cannot replace example input.")
        sys.exit(1)

    example_path = Path(config.predictkcat_DLKcat) / "DLKcat" / \
        "DeeplearningApproach" / "Code" / "example" / "input.tsv"
    example_parent = example_path.parent
    if not example_parent.exists():
        print(f"→ Creating directory for example input at {example_parent}")
        example_parent.mkdir(parents=True, exist_ok=True)

    print(f"→ Copying {src} to {example_path}")
    shutil.copy2(src, example_path)
    print("→ Replacement complete.")

def unzip_input_zip():
    zip_path = Path(config.predictkcat_DLKcat) / "DLKcat" / \
        "DeeplearningApproach" / "Data" / "input.zip"
    extract_dir = zip_path.parent
    if not zip_path.exists():
        print(f"→ Zip file {zip_path} not found; cannot unzip.")
        sys.exit(1)

    print(f"→ Unzipping {zip_path} to {extract_dir}")
    with zipfile.ZipFile(zip_path, 'r') as z:
        z.extractall(extract_dir)
    print("→ Unzip complete.")

def run_prediction():
    example_dir = Path(config.predictkcat_DLKcat) / "DLKcat" / \
        "DeeplearningApproach" / "Code" / "example"
    script     = example_dir / "prediction_for_input.py"
    input_file = example_dir / "input.tsv"

    if not script.exists():
        print(f"→ Prediction script {script} not found.")
        sys.exit(1)

    if not input_file.exists():
        print(f"→ Input file {input_file} not found.")
        sys.exit(1)

    print(f"→ Running prediction script in {example_dir}")
    # pass the input file path as argv[1]
    proc = subprocess.Popen(
        [sys.executable, str(script), str(input_file)],
        cwd=example_dir,
        stdout=sys.stdout,
        stderr=sys.stderr
    )
    ret = proc.wait()
    if ret != 0:
        print(f"→ Prediction script failed with exit code {ret}")
        sys.exit(ret)
    print("→ Prediction script finished successfully.")


def copy_output():
    src = Path(config.predictkcat_DLKcat) / "DLKcat" / \
        "DeeplearningApproach" / "Code" / "example" / "output.tsv"
    dst = Path(config.predictkcat_output_tsv)
    if not src.exists():
        print(f"→ Output file {src} not found; cannot copy.")
        sys.exit(1)

    dst_parent = dst.parent
    if not dst_parent.exists():
        dst_parent.mkdir(parents=True, exist_ok=True)

    print(f"→ Copying result from {src} to {dst}")
    shutil.copy2(src, dst)
    print("→ Output copy complete.")

def main():
    steps = [
        clone_dlkcat,
        replace_example_input,
        unzip_input_zip,
        run_prediction,
        copy_output,
    ]
    for step in steps:
        step()

if __name__ == "__main__":
    main()
