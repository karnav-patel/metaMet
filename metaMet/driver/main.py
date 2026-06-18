#!/usr/bin/env python3
import sys
import os
import time
import argparse
import subprocess
import tempfile
from pathlib import Path
import utils.installation as installation
import utils.precheck as precheck

def get_project_root():
    """
    Traverse upwards from the current file's location until the folder named 'metaMet' is found.
    """
    current_file_path = Path(__file__).resolve()
    for parent in current_file_path.parents:
        if parent.name == "metaMet":
            return parent
    raise RuntimeError("Project root folder 'metaMet' not found.")

def wait_for_job_completion(job_id, check_interval=60):
    """
    Polls the job status using squeue every `check_interval` seconds, printing progress updates.
    When the job is no longer listed, assumes completion.
    """
    progress = 0
    while True:
        result = subprocess.run(["squeue", "-j", str(job_id)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if str(job_id) not in result.stdout:
            print("Progress: 100%")
            break
        print(f"Progress: {progress}% (Job {job_id} still running)")
        time.sleep(check_interval)
        progress = min(progress + 10, 99)

def run_script_with_progress(script_path: Path):
    """
    Runs the given Python script as a subprocess while printing progress updates every minute.
    """
    print(f"Running {script_path} directly...")
    process = subprocess.Popen(["python", str(script_path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    progress = 0
    while process.poll() is None:
        print(f"Progress: {progress}% (Running {script_path.name})")
        time.sleep(60)
        progress = min(progress + 10, 99)
    # When finished, print any remaining output
    stdout, stderr = process.communicate()
    if stdout:
        print(stdout)
    if stderr:
        print(stderr)
    print(f"Progress: 100% (Finished {script_path.name})")

def submit_batch_job(scripts, job_name: str, job_dir: Path | None):
    """
    Creates one SLURM script that runs *all* Python files in order,
    then submits it via sbatch and waits for completion.
    """
    temp_dir = tempfile.gettempdir()
    sh_path = os.path.join(temp_dir, f"{job_name}_batch.sh")
    with open(sh_path, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"#SBATCH --job-name={job_name}\n")
        f.write("#SBATCH --partition=gpu\n")
        f.write("#SBATCH --gres=gpu:1\n")
        f.write("#SBATCH --mem=64G\n")
        f.write("#SBATCH --cpus-per-task=1\n")
        f.write("#SBATCH --time=100:00:00\n")
        f.write(f"#SBATCH --output=/home/yaolab/kpatel/download/job/logs/{job_name}_batch_output.log\n")
        f.write(f"#SBATCH --error=/home/yaolab/kpatel/download/job/logs/{job_name}_batch_error.log\n\n")
        f.write("module load cuda\n")
        f.write("module load anaconda\n")
        f.write("conda activate python\n\n")
        if job_dir is not None:
            f.write(f'cd "{job_dir}"\n\n')
        # add one `python …` line per script, preserving order
        for script in scripts:
            f.write(f'python "{script}"\n')
    os.chmod(sh_path, 0o755)

    print(f"Submitting batch job script: {sh_path}")
    result = subprocess.run(["sbatch", sh_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print(result.stdout)
    if result.stderr:
        print(result.stderr)
    # grab job ID and wait
    job_id = next((tok for tok in result.stdout.split() if tok.isdigit()), None)
    if job_id:
        print(f"Batch Job ID: {job_id}")
        wait_for_job_completion(job_id)
    else:
        print("Failed to retrieve job ID from sbatch output.")

def resolve_base_dir(project_root: Path, use_modeling: bool) -> Path:
    """
    Returns the selected base directory under the project root.
    """
    return project_root / ("data_modeling" if use_modeling else "data_preprocessing")

def resolve_scripts(base_dir: Path, folder: str | None, scripts_cli: list[str] | None, types: list[str] | None):
    """
    Build a list of script Paths to run, based on:
      - Explicit scripts via -script (highest priority; order preserved)
      - Or discovery via -type within the provided -folder
    """
    resolved = []

    # If explicit scripts are provided, honor them and ignore -type filters.
    if scripts_cli:
        # Determine the anchor directory for resolving relatives
        anchor_dir = base_dir / folder if folder else base_dir
        for s in scripts_cli:
            p = Path(s)
            if not p.is_absolute():
                # try under folder (if given) or base_dir
                p = (anchor_dir / s).resolve()
            resolved.append(p)
        return resolved

    # Otherwise, use type-based discovery (requires folder)
    if not folder:
        raise ValueError("When using -type filters, you must also provide -folder to search within.")

    job_dir = (base_dir / folder)
    if not job_dir.exists():
        raise FileNotFoundError(f"The folder {job_dir} does not exist.")

    if any(t.lower() == "all" for t in (types or [])):
        resolved = [f.resolve() for f in job_dir.glob("*.py")]
    else:
        resolved_set = set()
        for t in (types or []):
            for f in job_dir.glob(f"*{t}*.py"):
                resolved_set.add(f.resolve())
        resolved = sorted(resolved_set, key=lambda p: p.name)

    return resolved

def main():
    project_root = get_project_root()

    epilog = """
Examples:
  # Init only (installs + precheck using config in selected base; defaults to preprocessing)
  meta_runner.py -init

  # Run locally all *.py in a folder under data_preprocessing matching types 'ec' or 'reaction'
  meta_runner.py -run -p -folder job1 -type ec -type reaction

  # Submit a SLURM job for every *.py in a modeling folder
  meta_runner.py -job my_kcat_job -m -folder kcat_stage -type all

  # Run a single script locally (relative to base folder)
  meta_runner.py -run -p -script prepare_ecs.py

  # Submit a SLURM job running two explicit scripts (relative to a subfolder)
  meta_runner.py -job batch_A -m -folder stage2 -script step1_fit.py -script step2_eval.py

  # Submit a SLURM job using absolute script paths (no folder needed)
  meta_runner.py -job freeform -script /path/to/foo.py -script /path/to/bar.py
"""

    parser = argparse.ArgumentParser(
        description="Run data scripts in local or SLURM modes across data_preprocessing or data_modeling.",
        epilog=epilog,
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Base selection
    base_group = parser.add_mutually_exclusive_group()
    base_group.add_argument("-m", "--modeling", action="store_true",
                            help="Use 'data_modeling' as the base directory.")
    base_group.add_argument("-p", "--preprocessing", action="store_true",
                            help="Use 'data_preprocessing' as the base directory (default).")

    # Init
    parser.add_argument("-init", action="store_true",
                        help="Run installation and precheck steps using the config file in the selected base directory.")

    # Mode selection
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("-job", type=str,
                            help=("Submit as SLURM job with this JOB NAME. "
                                  "Use -folder with -type to discover scripts, or -script to specify exact files."))
    mode_group.add_argument("-run", action="store_true",
                            help="Run scripts directly in the terminal.")

    # What to run
    parser.add_argument("-folder", type=str,
                        help=("Optional subfolder under the selected base in which to search for scripts "
                              "(used with -type) or to resolve relative -script paths."))
    parser.add_argument("-type", action="append",
                        help=("Type filter(s) to select scripts within the -folder. "
                              "Use '-type all' to run all scripts in the folder; "
                              "otherwise provide substrings like '-type ec' or '-type reaction'."))
    parser.add_argument("-script", action="append",
                        help=("Path to a Python script to run (can repeat). "
                              "Relative paths are resolved under the selected base and optional -folder. "
                              "When provided, -type filters are ignored."))

    args = parser.parse_args()

    # Choose base; default to preprocessing if neither -m nor -p is given
    use_modeling = bool(args.modeling)
    base_dir = resolve_base_dir(project_root, use_modeling)

    # If -init is specified, run installation and precheck.
    if args.init:
        # Prefer config within selected base; fallback to preprocessing config if missing
        candidate_cfg = base_dir / "config" / "config.py"
        if not candidate_cfg.exists():
            candidate_cfg = project_root / "data_preprocessing" / "config" / "config.py"

        print(f"Using config file: {candidate_cfg}")
        print("Running installation...")
        installation.main()  # Installs requirements and creates necessary directories.
        print("Running precheck...")
        sys.argv = [sys.argv[0], "--config", str(candidate_cfg)]
        precheck.main()
        print("Precheck passed. Installation and precheck steps completed.")
        # If -init is the only action, exit now.
        if not args.job and not args.run and not args.type and not args.script and not args.folder:
            return

    # If not init-only, enforce mode selection and target specification
    if not args.job and not args.run:
        parser.error("Either -job or -run must be provided (unless using only -init).")

    # You must specify at least one of -script or -type
    if not args.script and not args.type:
        parser.error("Provide at least one of -script or -type (unless using only -init).")

    # If using -type without -script, require -folder to know where to search
    if args.type and not args.script and not args.folder:
        parser.error("When using -type filters, -folder is required to indicate the search directory.")

    # Determine job name and working dir anchor (for batch)
    if args.job:
        job_name = args.job
    else:
        # local run: invent a name for logging/prints only
        job_name = args.folder or "local_run"

    # Establish a job_dir anchor:
    # - If folder provided: base_dir/folder
    # - Else: base_dir (so relative -script paths still resolve cleanly)
    job_dir = (base_dir / args.folder) if args.folder else base_dir

    # Resolve which scripts to run
    try:
        script_files = resolve_scripts(base_dir, args.folder, args.script, args.type)
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Validate files exist
    missing = [str(p) for p in script_files if not p.exists()]
    if missing:
        print("The following script(s) do not exist:")
        for m in missing:
            print(f"  - {m}")
        sys.exit(1)

    # If discovery via -type (no -script), sort by filename for deterministic order.
    # If -script was provided, preserve CLI order.
    if not args.script:
        script_files = sorted(script_files, key=lambda p: p.name)

    print("Scripts to run:")
    for s in script_files:
        print(f" - {s}")

    # Execute
    if args.job:
        submit_batch_job(script_files, job_name, job_dir)
    else:
        for script in script_files:
            run_script_with_progress(script)

    print("All tasks completed.")

if __name__ == "__main__":
    main()
