import subprocess
import sys
import pandas as pd
from pathlib import Path


def test_build_funnel_whitelist(tmp_path):
    # Mock a Glide CSV
    g_csv = tmp_path / "glide_dock_f1.csv"
    g_csv.write_text(
        "Title,docking_score\nMOL1,-8.5\nMOL2,-7.0\nMOL3,-5.0\nCTRL_1,-9.0\n"
    )

    # Mock a GNINA SDF
    n_sdf = tmp_path / "gnina_score_only_f1.sdf"
    n_sdf.write_text(
        """MOL1
$$$$
> <CNNaffinity>
6.5

$$$$
MOL2
$$$$
> <CNNaffinity>
5.0

$$$$
MOL3
$$$$
> <CNNaffinity>
4.0

$$$$
"""
    )

    out_csv = tmp_path / "whitelist.csv"
    script = Path(__file__).parent.parent.parent / "bin" / "build_funnel_whitelist.py"

    cmd = [
        sys.executable,
        str(script),
        "--glide_csvs",
        str(g_csv),
        "--gnina_sdfs",
        str(n_sdf),
        "--top_k_frames",
        "1",
        "--percentile_cutoff",
        "0.50",  # Top 50%
        "--benchmark_ligands",
        "CTRL_1",
        "--output",
        str(out_csv),
    ]

    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"build_funnel failed: {res.stderr}"

    df = pd.read_csv(out_csv)
    titles = set(df["Title"].tolist())

    # We expect CTRL_1 (forced) and MOL1 (best triage score)
    assert "CTRL_1" in titles
    assert "MOL1" in titles
    # MOL3 should be dropped (bottom 50%)
    assert "MOL3" not in titles
