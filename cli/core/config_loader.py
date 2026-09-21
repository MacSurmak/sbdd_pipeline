"""
File: cli/core/config_loader.py
Description: Scans, loads, and merges atomic configuration YAML files from configs/.
"""
from pathlib import Path
from typing import Dict, Any
import yaml

def load_project_configs(project_dir: Path) -> Dict[str, Any]:
    """Merges all YAML files in project_dir/configs or project_dir into a single dictionary."""
    config_dir = project_dir / "configs"
    search_dirs = [config_dir] if config_dir.exists() else [project_dir]
    
    merged: Dict[str, Any] = {}
    yaml_files = []
    for sdir in search_dirs:
        yaml_files.extend(list(sdir.glob("*.yaml")) + list(sdir.glob("*.yml")))

    for yf in sorted(yaml_files):
        try:
            with open(yf, "r", encoding="utf-8") as f:
                data = yaml.safe_load(f)
                if isinstance(data, dict):
                    merged.update(data)
        except Exception as e:
            print(f"Warning: Could not parse {yf.name}: {e}")

    # Auto-detection fallbacks for inputs in project root
    if "receptor" not in merged or not merged["receptor"].get("input_file"):
        cif_files = list(project_dir.glob("*.cif"))
        pdb_files = list(project_dir.glob("*.pdb"))
        if cif_files:
            merged.setdefault("receptor", {})["input_file"] = cif_files[0].name
        elif pdb_files:
            merged.setdefault("receptor", {})["input_file"] = pdb_files[0].name

    if "screening" not in merged or not merged["screening"].get("library_file"):
        sdf_files = [f for f in project_dir.glob("*.sdf") if "dock" not in f.name.lower()]
        if sdf_files:
            merged.setdefault("screening", {})["library_file"] = sdf_files[0].name

    return merged

def build_nextflow_args(config: Dict[str, Any], project_dir: Path) -> list:
    """Translates dictionary of parameters into Nextflow CLI arguments."""
    args = []

    # Project
    proj_cfg = config.get("project", {})
    proj_name = proj_cfg.get("name", project_dir.name)
    args.extend(["--project_name", proj_name])

    # Receptor
    rec_cfg = config.get("receptor", {})
    if rec_cfg.get("input_file"):
        rec_path = project_dir / rec_cfg["input_file"]
        if not rec_path.exists():
            # Check receptor/ folder
            rec_path = project_dir / "receptor" / rec_cfg["input_file"]
        args.extend(["--receptor", str(rec_path.resolve())])

    if rec_cfg.get("native_ligand_resname"):
        args.extend(["--native_ligand", str(rec_cfg["native_ligand_resname"])])

    if rec_cfg.get("model_index") is not None:
        args.extend(["--model_index", str(rec_cfg["model_index"])])

    if rec_cfg.get("keep_cofactors"):
        cofactors = rec_cfg["keep_cofactors"]
        if isinstance(cofactors, list):
            args.extend(["--keep_cofactors", " ".join(cofactors)])

    # Membrane
    mem_cfg = config.get("membrane", {})
    if "enabled" in mem_cfg:
        args.extend(["--membrane", "true" if mem_cfg["enabled"] else "false"])
    if mem_cfg.get("type"):
        args.extend(["--membrane_type", str(mem_cfg["type"])])
    if mem_cfg.get("topology"):
        args.extend(["--membrane_topo", str(mem_cfg["topology"])])
    if mem_cfg.get("chains"):
        chains = mem_cfg["chains"]
        args.extend(["--membrane_chains", ",".join(chains) if isinstance(chains, list) else str(chains)])

    # Library
    screen_cfg = config.get("screening", {})
    if screen_cfg.get("library_file"):
        lib_path = project_dir / screen_cfg["library_file"]
        if not lib_path.exists():
            lib_path = project_dir / "ligands" / screen_cfg["library_file"]
        if lib_path.exists():
            args.extend(["--library", str(lib_path.resolve())])

    # Docking Engine
    gnina_cfg = config.get("gnina", {})
    glide_cfg = config.get("glide", {})
    if glide_cfg.get("enabled", False):
        args.extend(["--dock_engine", "glide"])
    else:
        args.extend(["--dock_engine", "gnina"])
        if gnina_cfg.get("exhaustiveness"):
            args.extend(["--exhaustiveness", str(gnina_cfg["exhaustiveness"])])
        if gnina_cfg.get("cnn_scoring"):
            args.extend(["--cnn_scoring", str(gnina_cfg["cnn_scoring"])])
        if gnina_cfg.get("batch_size"):
            args.extend(["--batch_size", str(gnina_cfg["batch_size"])])

    # MD Settings
    md_cfg = config.get("md_simulation", {})
    if "enabled" in md_cfg:
        args.extend(["--run_md", "true" if md_cfg["enabled"] else "false"])
    prod_cfg = md_cfg.get("production", {})
    if prod_cfg.get("total_time_ns"):
        args.extend(["--md_total_ns", str(prod_cfg["total_time_ns"])])
    if prod_cfg.get("chunk_size_ns"):
        args.extend(["--md_chunk_ns", str(prod_cfg["chunk_size_ns"])])
    if prod_cfg.get("replicates"):
        args.extend(["--md_replicates", str(prod_cfg["replicates"])])
    if prod_cfg.get("temperature"):
        args.extend(["--md_temperature", str(prod_cfg["temperature"])])

    # Output directory
    args.extend(["--outdir", str(project_dir.resolve())])

    return args

