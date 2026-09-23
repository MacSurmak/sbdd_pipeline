import subprocess
import sys
import csv
from pathlib import Path


def test_merge_csv(tmp_path):
    csv1 = tmp_path / "1.csv"
    csv1.write_text("Title,Score\nMOL1,1.5\nMOL2,2.0\n")

    csv2 = tmp_path / "2.csv"
    csv2.write_text("Title,Score\nMOL3,3.5\n")

    out_csv = tmp_path / "merged.csv"
    script = Path(__file__).parent.parent.parent / "bin" / "merge_csv.py"

    cmd = [
        sys.executable,
        str(script),
        "--inputs",
        str(csv1),
        str(csv2),
        "--output",
        str(out_csv),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"merge_csv failed: {res.stderr}"

    with open(out_csv, "r") as f:
        reader = list(csv.reader(f))

    assert len(reader) == 4  # 1 header + 3 data rows
    assert reader[0] == ["Title", "Score"]
    assert ["MOL1", "1.5"] in reader
    assert ["MOL3", "3.5"] in reader
