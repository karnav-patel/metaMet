#!/usr/bin/env python3
import sys
import subprocess
from pathlib import Path

def get_project_root():
    """
    Traverse upwards from the current file's location until the folder named 'metaMet' is found.
    This ensures that even if this file is in a subfolder (e.g., utils/), we locate the actual project root.
    """
    current_file_path = Path(__file__).resolve()
    for parent in current_file_path.parents:
        if parent.name == "metaMet":
            return parent
    raise RuntimeError("Project root folder 'metaMet' not found.")

def install_requirements(requirements_file: Path) -> bool:
    if not requirements_file.exists():
        print(f"[ERROR] Requirements file not found: {requirements_file}")
        return False
    print(f"[INFO] Installing Python packages from: {requirements_file}")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", str(requirements_file)])
        print("[OK] Python packages installed successfully.")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to install Python packages. Error: {e}")
        return False

def create_directory(dir_path: Path, description: str):
    if dir_path.exists():
        print(f"[OK] {description} already exists: {dir_path}")
    else:
        try:
            dir_path.mkdir(parents=True, exist_ok=True)
            print(f"[OK] Created {description}: {dir_path}")
        except Exception as e:
            print(f"[ERROR] Failed to create {description} at {dir_path}. Error: {e}")

def main():
    # Determine the project root using the helper function.
    project_root = get_project_root()
    requirements_file = project_root / "driver" / "requirements.txt"

    # Install required Python packages.
    if not install_requirements(requirements_file):
        sys.exit(1)

    # List of directories to create with a description.
    directories = [
        ("BRenda raw data directory", project_root / "data" / "raw" / "brenda"),
        ("MetaCyc raw data directory", project_root / "data" / "raw" / "metacyc"),
        ("MetaCyc raw data directory", project_root / "data" / "raw" / "mapping_ec_numbers"),
        ("BRenda processed data directory", project_root / "data" / "processed" / "brenda"),
        ("KEGG processed data directory", project_root / "data" / "processed" / "kegg"),
        ("KEGG more processed data directory", project_root / "data" / "processed" / "kegg" / "more"),
        ("MetaCyc processed data directory", project_root / "data" / "processed" / "metacyc"),
        ("MetaCyc more processed data directory", project_root / "data" / "processed" / "metacyc" / "more"),
        ("Merged processed data directory", project_root / "data" / "processed" / "merged"),
        ("Sabio processed data directory", project_root / "data" / "processed" / "sabio"),
        ("Sabio more processed data directory", project_root / "data" / "processed" / "sabio" / "more"),
        ("PredictKcat processed data directory", project_root / "data" / "processed" / "predictkcat"),
        ("PredictKcat more processed data directory", project_root / "data" / "processed" / "predictkcat" / "more"),
        ("overview processed data directory", project_root / "data" / "processed" / "overview"),
    ]

    for desc, path in directories:
        create_directory(path, desc)

    print("\nInstallation completed successfully.")

if __name__ == "__main__":
    main()
