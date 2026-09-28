import subprocess
import sys
from pathlib import Path


def test_select_hits(tmp_path):
    csv_file = tmp_path / "ranking.csv"
    csv_file.write_text(
        "Rank,Title,GlideScore,Dock_CNNaffinity\n1,PPZ_1H_Ph_BnMe,-8.5,6.2\n2,PPZ_1H_Ph_Pr2,-7.8,5.9\n"
    )
    out_file = tmp_path / "hits.txt"

    script = Path(__file__).parent.parent.parent / "bin" / "select_hits.py"
    cmd = [
        sys.executable,
        str(script),
        "--ranking_csv",
        str(csv_file),
        "--top_n",
        "1",
        "--include_native",
        "true",
        "--native_ligand",
        "BTN",
        "--out",
        str(out_file),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0
    lines = [l.strip() for l in out_file.read_text().splitlines() if l.strip()]
    assert lines == ["PPZ_1H_Ph_BnMe,1", "BTN,1"]
