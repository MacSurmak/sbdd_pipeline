import subprocess
import sys
from pathlib import Path

SAMPLE_GRO = """Declash Test System
    4
    1MOL1     C1    1   1.000   1.000   1.000
    1MOL1     C2    2   1.100   1.000   1.000
    2MOL2     C1    3   1.020   1.000   1.000
    2MOL2     C2    4   1.120   1.000   1.000
   5.00000   5.00000   5.00000
"""


def make_gro_line(resid, resname, atomname, atomnr, x, y, z):
    return f"{resid:5d}{resname[:5]:<5}{atomname[:5]:>5}{atomnr:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n"


SAMPLE_TOP = """[ moleculetype ]
; Name   nrexcl
MOL1     3

[ atoms ]
; nr  type  resnr resid atom cgnr charge mass
1     c3    1     MOL1  C1   1    0      12.011
2     c3    1     MOL1  C2   1    0      12.011

[ moleculetype ]
; Name   nrexcl
MOL2     3

[ atoms ]
; nr  type  resnr resid atom cgnr charge mass
1     c3    2     MOL2  C1   1    0      12.011
2     c3    2     MOL2  C2   1    0      12.011

[ system ]
Test

[ molecules ]
MOL1  1
MOL2  1
"""


def test_declash_gmx(tmp_path):
    gro_content = "Declash Test System\n    4\n"
    gro_content += make_gro_line(1, "MOL1", "C1", 1, 1.000, 1.000, 1.000)
    gro_content += make_gro_line(1, "MOL1", "C2", 2, 1.100, 1.000, 1.000)
    gro_content += make_gro_line(2, "MOL2", "C1", 3, 1.020, 1.000, 1.000)
    gro_content += make_gro_line(2, "MOL2", "C2", 4, 1.120, 1.000, 1.000)
    gro_content += "   5.00000   5.00000   5.00000\n"

    gro_file = tmp_path / "clashed.gro"
    gro_file.write_text(gro_content)

    top_file = tmp_path / "topol.top"
    top_file.write_text(SAMPLE_TOP)

    out_gro = tmp_path / "declashed.gro"
    report_txt = tmp_path / "report.txt"

    declash_script = Path(__file__).parent.parent.parent / "bin" / "declash_gmx.py"

    cmd = [
        sys.executable,
        str(declash_script),
        "-f",
        str(gro_file),
        "-p",
        str(top_file),
        "-o",
        str(out_gro),
        "--cutoff",
        "0.15",
        "--niter",
        "50",
        "--report",
        str(report_txt),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"declash failed: {res.stderr}\n{res.stdout}"

    assert out_gro.exists()
    assert report_txt.exists()

    # Read declashed gro and verify coordinates moved apart
    lines = out_gro.read_text().splitlines()
    assert len(lines) == 7
    # Original distance was 0.02 nm between atom 1 and 3 (1.000 vs 1.020).
    # After declash with cutoff 0.15 nm, distance should increase.
    x1 = float(lines[2][20:28])
    x3 = float(lines[4][20:28])
    assert abs(x3 - x1) > 0.03
