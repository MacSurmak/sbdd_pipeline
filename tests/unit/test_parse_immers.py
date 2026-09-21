import json
import subprocess
import sys
from pathlib import Path

def test_parse_immers(sample_datapar1, tmp_path):
    out_json = tmp_path / "membrane.json"
    out_txt = tmp_path / "membrane.txt"
    parse_script = Path(__file__).parent.parent.parent / "bin" / "parse_immers.py"

    cmd = [
        sys.executable,
        str(parse_script),
        "--datapar", str(sample_datapar1),
        "--out_json", str(out_json),
        "--out_txt", str(out_txt)
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"parse_immers failed: {res.stderr}"

    assert out_json.exists()
    assert out_txt.exists()

    with open(out_json) as f:
        data = json.load(f)

    assert data["thickness"] == 20.7
    assert data["tilt"] == 28.0
    assert data["dg"] == -20.2
    assert out_txt.read_text().strip() == "20.7"

