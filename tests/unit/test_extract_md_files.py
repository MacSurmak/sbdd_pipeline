import subprocess
import sys
from pathlib import Path

SAMPLE_PDB = """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.362   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.272   2.348   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.966  -0.764  -1.217  1.00  0.00           C
TER
END
"""


def test_extract_md_files_rdkit(sample_sdf, tmp_path):
    rec_pdb = tmp_path / "receptor.pdb"
    rec_pdb.write_text(SAMPLE_PDB)

    outdir = tmp_path / "extracted"
    extract_script = Path(__file__).parent.parent.parent / "bin" / "extract_md_files.py"

    cmd = [
        sys.executable,
        str(extract_script),
        "--lig_file",
        str(sample_sdf),
        "--lig_title",
        "MOL1",
        "--rec_file",
        str(rec_pdb),
        "--out_dir",
        str(outdir),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"extract_md_files failed: {res.stderr}\n{res.stdout}"

    assert (outdir / "receptor.pdb").exists()
    assert (outdir / "lig_charge.txt").exists()
    assert (outdir / "ligand.mol2").exists() or (outdir / "ligand.sdf").exists()
