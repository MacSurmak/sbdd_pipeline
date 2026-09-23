import json
import subprocess
import sys
from pathlib import Path
import pandas as pd


def test_summarize_results_gnina(sample_sdf, tmp_path):
    docking_dir = tmp_path / "docking"
    docking_dir.mkdir()
    # Put sample_sdf into docking_dir matching canonical gnina dock naming
    (docking_dir / "test_f1_gnina_dock.sdf").write_text(sample_sdf.read_text())

    reports_dir = tmp_path / "reports"
    summarize_script = (
        Path(__file__).parent.parent.parent / "bin" / "summarize_results.py"
    )

    cmd = [
        sys.executable,
        str(summarize_script),
        "--docking_dir",
        str(docking_dir),
        "--outdir",
        str(reports_dir),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"summarize_results failed: {res.stderr}"

    ranking_csv = reports_dir / "ranking.csv"
    summary_json = reports_dir / "summary.json"

    assert ranking_csv.exists()
    assert summary_json.exists()

    df = pd.read_csv(ranking_csv)
    assert len(df) == 3
    # Sorted by Dock_CNNaffinity descending: MOL3 (8.4) > MOL1 (7.2) > MOL2 (6.1)
    assert df.iloc[0]["Title"] == "MOL3"
    assert df.iloc[0]["Dock_CNNaffinity"] == 8.4
    assert df.iloc[1]["Title"] == "MOL1"
    assert df.iloc[2]["Title"] == "MOL2"
