import subprocess
import sys
from pathlib import Path

SAMPLE_PDB_RAW = """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.362   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.272   2.348   0.000  1.00  0.00           O
ATOM      5  H1  ALA A   1      -0.500   0.000   0.800  1.00  0.00           H
ATOM      6  N   NMA A   2       3.200   1.362   0.000  1.00  0.00           N
ATOM      7  H   NMA A   2       3.700   0.500   0.000  1.00  0.00           H
ATOM      8  CH3 NMA A   2       3.900   2.500   0.000  1.00  0.00           C
TER
END
"""

SAMPLE_MOL2 = """@<TRIPOS>MOLECULE
UNK
 2 1 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1          5.0000    5.0000    5.0000 c3        1 UNK       0.0000
      2 C2          6.2000    5.0000    5.0000 c3        1 UNK       0.0000
@<TRIPOS>BOND
     1     1     2 1
"""

def test_amber_fix_and_merge(tmp_path):
    raw_pdb = tmp_path / "raw.pdb"
    raw_pdb.write_text(SAMPLE_PDB_RAW)

    fixed_pdb = tmp_path / "fixed.pdb"
    lig_mol2 = tmp_path / "lig.mol2"
    lig_mol2.write_text(SAMPLE_MOL2)
    complex_pdb = tmp_path / "complex.pdb"

    amber_fix_script = Path(__file__).parent.parent.parent / "bin" / "amber_fix_pdb.py"
    merge_script = Path(__file__).parent.parent.parent / "bin" / "merge_complex.py"

    # 1. amber_fix_pdb
    cmd_fix = [sys.executable, str(amber_fix_script), str(raw_pdb), str(fixed_pdb)]
    res_fix = subprocess.run(cmd_fix, capture_output=True, text=True)
    assert res_fix.returncode == 0, f"amber_fix failed: {res_fix.stderr}"
    assert fixed_pdb.exists()

    fixed_text = fixed_pdb.read_text()
    # Check H1 renamed to H, NMA renamed to NME
    assert "NME" in fixed_text

    # 2. merge_complex
    cmd_merge = [sys.executable, str(merge_script), str(fixed_pdb), str(lig_mol2), str(complex_pdb)]
    res_merge = subprocess.run(cmd_merge, capture_output=True, text=True)
    assert res_merge.returncode == 0, f"merge_complex failed: {res_merge.stderr}"
    assert complex_pdb.exists()

    complex_text = complex_pdb.read_text()
    assert "UNK" in complex_text
    assert "ALA" in complex_text

