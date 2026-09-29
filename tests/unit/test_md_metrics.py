import json
import subprocess
import sys
from pathlib import Path

SAMPLE_PDB_COMPLEX = """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.362   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.272   2.348   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.966  -0.764  -1.217  1.00  0.00           C
TER
HETATM    6  C1  UNK Z   1       3.000   2.000   0.000  1.00  0.00           C
HETATM    7  C2  UNK Z   1       4.200   2.000   0.000  1.00  0.00           C
TER
END
"""


def test_md_metrics_calculation(tmp_path):
    pdb_file = tmp_path / "complex.pdb"
    pdb_file.write_text(SAMPLE_PDB_COMPLEX)

    outdir = tmp_path / "metrics"
    metrics_script = Path(__file__).parent.parent.parent / "bin" / "md_metrics.py"

    cmd = [
        sys.executable,
        str(metrics_script),
        "--pdb",
        str(pdb_file),
        "--xtc",
        str(pdb_file),
        "--outdir",
        str(outdir),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"md_metrics failed: {res.stderr}\n{res.stdout}"

    assert (outdir / "rmsd.csv").exists()
    assert (outdir / "rmsf.csv").exists()
    assert (outdir / "contacts.csv").exists()
    assert (outdir / "md_metrics_summary.json").exists()

    with open(outdir / "md_metrics_summary.json") as f:
        summary = json.load(f)

    assert "Prot_BB_RMSD_Mean_A" in summary
    assert "Lig_Heavy_RMSD_Mean_A" in summary
