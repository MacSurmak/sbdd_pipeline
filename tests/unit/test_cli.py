import subprocess
import sys
from pathlib import Path
from cli.core.config_loader import load_project_configs, build_nextflow_args


def test_config_loader(tmp_path):
    configs_dir = tmp_path / "configs"
    configs_dir.mkdir()

    (tmp_path / "test_rec.pdb").write_text("ATOM 1 ...\n")
    (configs_dir / "project.yaml").write_text("project:\n  name: TestProj\n")
    (configs_dir / "receptor.yaml").write_text(
        "receptor:\n  input_file: test_rec.pdb\n  native_ligand_resname: PKA\n"
    )
    (configs_dir / "gnina.yaml").write_text(
        "gnina:\n  enabled: true\n  exhaustiveness: 8\n  cnn_scoring: rescore\n  batch_size: 10\n"
    )
    (configs_dir / "membrane.yaml").write_text("membrane:\n  enabled: false\n")

    cfg = load_project_configs(tmp_path)
    assert cfg["project"]["name"] == "TestProj"
    assert cfg["gnina"]["exhaustiveness"] == 8

    nf_args = build_nextflow_args(cfg, tmp_path)
    assert "--project_name" in nf_args
    assert "TestProj" in nf_args
    assert "--exhaustiveness" in nf_args
    assert "8" in nf_args
    assert "--membrane" in nf_args
    assert "false" in nf_args


def test_sbdd_init(tmp_path):
    import shutil

    target = tmp_path / "new_project"
    sbdd_bin = shutil.which("sbdd") or str(Path(sys.prefix) / "bin" / "sbdd")
    cmd = [sbdd_bin, "init", str(target)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"sbdd init failed: {res.stderr}"

    assert (target / "configs").exists()
    assert (target / "configs" / "project.yaml").exists()
    assert (target / "configs" / "gnina.yaml").exists()
    assert (target / "configs" / "receptor.yaml").exists()
    assert (target / "configs" / "membrane.yaml").exists()
