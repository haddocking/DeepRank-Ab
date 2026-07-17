"""
Smoke test for the end-to-end inference pipeline (scripts/inference.py).

Runs the pipeline against example/test.pdb and checks that a prediction CSV
with the expected columns is produced. Meant to run inside the Docker image
(see Dockerfile) where all heavy dependencies (ANARCI, ESM, torch, pdb2sql,
voronota, freesasa...) are installed.
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
TEST_PDB = REPO_ROOT / "example" / "test.pdb"


def test_inference_smoke(tmp_path):
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts" / "inference.py"), str(TEST_PDB)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=1800,
    )
    # pipeline ran to completion without error
    assert result.returncode == 0, result.stdout + result.stderr

    # workspace dir for the run was created (name encodes detected chain mapping)
    workdirs = list(tmp_path.glob("test-deeprank_ab_pred_*"))
    assert workdirs, (
        f"no workspace created.\nstdout={result.stdout}\nstderr={result.stderr}"
    )

    # final predictions CSV was written into that workspace
    csvs = list(workdirs[0].glob("*_predictions.csv"))
    assert csvs, f"no prediction CSV found in {workdirs[0]}"

    # CSV has the columns downstream consumers rely on
    df = pd.read_csv(csvs[0])
    for col in ("pdb_id", "predicted_dockq", "HL_contact_flag", "vdw_clash_flag"):
        assert col in df.columns
    # and at least one scored model
    assert len(df) >= 1
