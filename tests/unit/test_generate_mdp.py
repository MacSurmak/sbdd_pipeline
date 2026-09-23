import subprocess
import sys
from pathlib import Path


def test_generate_mdp(tmp_path):
    tpl_dir = tmp_path / "templates"
    tpl_dir.mkdir()

    (tpl_dir / "minim.mdp").write_text("integrator=steep\n{DEFINE}\n")
    (tpl_dir / "nvt.mdp").write_text("ref_t={TEMP}\nnsteps={NSTEPS}\n")
    (tpl_dir / "npt.mdp").write_text("pcoupltype={PCOUPL_TYPE}\nnsteps={NSTEPS}\n")

    out_dir = tmp_path / "out"
    out_dir.mkdir()

    script = Path(__file__).parent.parent.parent / "bin" / "generate_mdp.py"

    cmd = [
        sys.executable,
        str(script),
        "--templates_dir",
        str(tpl_dir),
        "--outdir",
        str(out_dir),
        "--is_membrane",
        "true",
        "--temperature",
        "310.15",
        "--nvt_short_ps",
        "50",
        "--nvt_long_ps",
        "500",
        "--npt_berendsen_ps",
        "1000",
        "--npt_pr_ps",
        "1000",
    ]

    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"generate_mdp failed: {res.stderr}"

    assert (out_dir / "stages_list.txt").exists()
    assert (out_dir / "nvt_short.mdp").exists()

    nvt_content = (out_dir / "nvt_short.mdp").read_text()
    assert "ref_t=310.15" in nvt_content
    assert "nsteps=50000" in nvt_content  # 50ps / 0.001dt = 50000

    npt_content = (out_dir / "npt_berendsen.mdp").read_text()
    assert "pcoupltype=semiisotropic" in npt_content
