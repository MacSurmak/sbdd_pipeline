import subprocess
import sys
from pathlib import Path

SAMPLE_NDX = """[ System ]
 1 2 3 4 5 6 7 8 9
[ Protein ]
 1 2 3 4
[ UNK ]
 5 6 7
[ SOL ]
 8 9
"""


def test_gromacs_ndx_fixer_ensure_complex(tmp_path):
    ndx_file = tmp_path / "index.ndx"
    ndx_file.write_text(SAMPLE_NDX)

    script = Path(__file__).parent.parent.parent / "bin" / "gromacs_ndx_fixer.py"

    # 1. Test ensure_complex
    cmd = [
        sys.executable,
        str(script),
        "--ndx",
        str(ndx_file),
        "--action",
        "ensure_complex",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0

    content = ndx_file.read_text()
    assert "[ COMPLEX ]" in content
    assert "1     2     3     4     5     6     7" in content


def test_gromacs_ndx_fixer_add_groups(tmp_path):
    ndx_file = tmp_path / "index.ndx"
    ndx_file.write_text(SAMPLE_NDX)

    script = Path(__file__).parent.parent.parent / "bin" / "gromacs_ndx_fixer.py"

    # 2. Test add_groups (Thermostat PROT_LIG / SOL_ION)
    cmd = [
        sys.executable,
        str(script),
        "--ndx",
        str(ndx_file),
        "--action",
        "add_groups",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0

    content = ndx_file.read_text()
    assert "[ PROT_LIG ]" in content
    assert "[ SOL_ION ]" in content

    # 8 and 9 are SOL, so SOL_ION should contain 8 and 9
    assert "8     9" in content
