import subprocess
import sys
from pathlib import Path

def test_split_sdf_basic(sample_sdf, tmp_path):
    outdir = tmp_path / "chunks"
    split_script = Path(__file__).parent.parent.parent / "bin" / "split_sdf.py"

    cmd = [
        sys.executable,
        str(split_script),
        "--input", str(sample_sdf),
        "--chunk_size", "2",
        "--outdir", str(outdir)
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"split_sdf failed: {res.stderr}"

    chunks = sorted(list(outdir.glob("chunk_*.sdf")))
    assert len(chunks) == 2, f"Expected 2 chunks for 3 molecules with batch_size=2, got {len(chunks)}"

    # Check chunk 0 has 2 molecules
    content_0 = chunks[0].read_text()
    assert content_0.count("$$$$") == 2
    assert "MOL1" in content_0
    assert "MOL2" in content_0

    # Check chunk 1 has 1 molecule
    content_1 = chunks[1].read_text()
    assert content_1.count("$$$$") == 1
    assert "MOL3" in content_1

