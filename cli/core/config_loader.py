"""
File: cli/core/config_loader.py
Description: Scans, loads, and strictly validates atomic configuration YAML files from configs/.
Zero silent fallbacks: required parameters must be explicitly specified.
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
            raise ValueError(f"Failed to parse YAML file {yf.name}: {e}")

    # Fallback to single inputs in project root if not explicitly in YAML
    if "receptor" not in merged or not merged["receptor"].get("input_file"):
        cif_files = list(project_dir.glob("*.cif"))
        pdb_files = list(project_dir.glob("*.pdb"))
        if cif_files:
            merged.setdefault("receptor", {})["input_file"] = cif_files[0].name
        elif pdb_files:
            merged.setdefault("receptor", {})["input_file"] = pdb_files[0].name

    if "screening" not in merged or not merged["screening"].get("library_file"):
        sdf_files = [
            f for f in project_dir.glob("*.sdf") if "dock" not in f.name.lower()
        ]
        if sdf_files:
            merged.setdefault("screening", {})["library_file"] = sdf_files[0].name

    return merged


def build_nextflow_args(config: Dict[str, Any], project_dir: Path) -> list:
    """Translates dictionary of parameters into Nextflow CLI arguments with strict validation."""
    args = []

    # 1. Project
    proj_cfg = config.get("project", {})
    if not proj_cfg.get("name"):
        raise ValueError(
            "Config validation failed: 'project.name' in project.yaml is required!"
        )
    args.extend(["--project_name", str(proj_cfg["name"])])

    # 2. Receptor
    rec_cfg = config.get("receptor", {})
    if not rec_cfg:
        raise ValueError(
            "Config validation failed: 'receptor' section in receptor.yaml is required!"
        )

    rec_file = rec_cfg.get("input_file")
    if not rec_file:
        raise ValueError(
            "Config validation failed: 'receptor.input_file' in receptor.yaml is required!"
        )
    rec_path = project_dir / rec_file
    if not rec_path.exists():
        rec_path = project_dir / "receptor" / rec_file
    if not rec_path.exists():
        raise FileNotFoundError(f"Receptor file not found: {rec_file} in {project_dir}")
    args.extend(["--receptor", str(rec_path.resolve())])

    if not rec_cfg.get("native_ligand_resname"):
        raise ValueError(
            "Config validation failed: 'receptor.native_ligand_resname' in receptor.yaml is required!"
        )
    args.extend(["--native_ligand", str(rec_cfg["native_ligand_resname"])])

    if "model_index" in rec_cfg:
        args.extend(["--model_index", str(rec_cfg["model_index"])])
    if "total_models" in rec_cfg:
        args.extend(["--total_models", str(rec_cfg["total_models"])])
    if rec_cfg.get("keep_cofactors"):
        cofactors = rec_cfg["keep_cofactors"]
        args.extend(
            [
                "--keep_cofactors",
                " ".join(cofactors) if isinstance(cofactors, list) else str(cofactors),
            ]
        )

    pw_cfg = rec_cfg.get("prep_wizard", {})
    if pw_cfg:
        if "cap_termini" in pw_cfg:
            args.extend(["--cap_termini", "true" if pw_cfg["cap_termini"] else "false"])
        if "target_ph" in pw_cfg:
            args.extend(["--target_ph", str(pw_cfg["target_ph"])])
        if "ph_tolerance" in pw_cfg:
            args.extend(["--ph_tolerance", str(pw_cfg["ph_tolerance"])])
        if "minimize_rmsd" in pw_cfg:
            args.extend(["--minimize_rmsd", str(pw_cfg["minimize_rmsd"])])

    # 3. Membrane
    mem_cfg = config.get("membrane", {})
    if not mem_cfg or "enabled" not in mem_cfg:
        raise ValueError(
            "Config validation failed: 'membrane.enabled' in membrane.yaml is required!"
        )

    if mem_cfg["enabled"]:
        args.extend(["--membrane", "true"])
        for key, flag in [("type", "--membrane_type"), ("topology", "--membrane_topo")]:
            if key not in mem_cfg:
                raise ValueError(
                    f"Config validation failed: 'membrane.{key}' is required when membrane is enabled!"
                )
            args.extend([flag, str(mem_cfg[key])])
        if "chains" in mem_cfg:
            chains = mem_cfg["chains"]
            args.extend(
                [
                    "--membrane_chains",
                    ",".join(chains) if isinstance(chains, list) else str(chains),
                ]
            )

        lipids_cfg = mem_cfg.get("lipids", {})
        if lipids_cfg:
            upper = lipids_cfg.get("upper_leaflet", {})
            lower = lipids_cfg.get("lower_leaflet", {})
            if "composition" in upper:
                args.extend(["--upper_lipids", str(upper["composition"])])
            if "ratio" in upper:
                args.extend(["--upper_ratio", str(upper["ratio"])])
            if "composition" in lower:
                args.extend(["--lower_lipids", str(lower["composition"])])
            if "ratio" in lower:
                args.extend(["--lower_ratio", str(lower["ratio"])])

        pad_cfg = mem_cfg.get("padding", {})
        if pad_cfg:
            if "lipid_padding" in pad_cfg:
                args.extend(["--lipid_padding", str(pad_cfg["lipid_padding"])])
            if "box_padding" in pad_cfg:
                args.extend(["--box_padding", str(pad_cfg["box_padding"])])

        ions_cfg = mem_cfg.get("ions", {})
        if ions_cfg and "salt_conc" in ions_cfg:
            args.extend(["--md_salt_conc", str(ions_cfg["salt_conc"])])
    else:
        args.extend(["--membrane", "false"])

    # 4. Library / Screening
    screen_cfg = config.get("screening", {})
    if screen_cfg.get("library_file"):
        lib_file = screen_cfg["library_file"]
        lib_path = project_dir / lib_file
        if not lib_path.exists():
            lib_path = project_dir / "ligands" / lib_file
        if not lib_path.exists():
            raise FileNotFoundError(
                f"Ligand library file not found: {lib_file} in {project_dir}"
            )
        args.extend(["--library", str(lib_path.resolve())])

    funnel_cfg = screen_cfg.get("funnel", {}) or config.get("analysis", {}).get(
        "funnel", {}
    )
    if funnel_cfg:
        if "top_frames_per_ligand" in funnel_cfg:
            args.extend(["--top_k_frames", str(funnel_cfg["top_frames_per_ligand"])])
        if "library_percentile_cutoff" in funnel_cfg:
            args.extend(
                ["--percentile_cutoff", str(funnel_cfg["library_percentile_cutoff"])]
            )
        if "benchmark_ligands" in funnel_cfg:
            benches = funnel_cfg["benchmark_ligands"]
            args.extend(
                [
                    "--benchmark_ligands",
                    " ".join(benches) if isinstance(benches, list) else str(benches),
                ]
            )

    # 5. Docking Engines (Glide / GNINA)
    gnina_cfg = config.get("gnina", {})
    glide_cfg = config.get("glide", {})
    if glide_cfg.get("enabled", False) and gnina_cfg.get("enabled", False):
        args.extend(["--dock_engine", "both"])
    elif glide_cfg.get("enabled", False):
        args.extend(["--dock_engine", "glide"])
    elif gnina_cfg.get("enabled", False):
        args.extend(["--dock_engine", "gnina"])

    if gnina_cfg.get("enabled", False):
        for key, flag in [
            ("exhaustiveness", "--exhaustiveness"),
            ("cnn_scoring", "--cnn_scoring"),
            ("num_modes", "--num_modes"),
            ("autobox_add", "--autobox_add"),
            ("cnn_model", "--cnn_model"),
            ("batch_size", "--batch_size"),
        ]:
            if key in gnina_cfg:
                args.extend([flag, str(gnina_cfg[key])])

    if glide_cfg.get("enabled", False):
        if "precision" in glide_cfg:
            args.extend(["--glide_precision", str(glide_cfg["precision"])])
        if "poses_per_lig" in glide_cfg:
            args.extend(["--poses_per_lig", str(glide_cfg["poses_per_lig"])])
        if "forcefield" in glide_cfg:
            args.extend(["--glide_forcefield", str(glide_cfg["forcefield"])])
        grid_cfg = glide_cfg.get("grid", {})
        if grid_cfg:
            if "inner_box" in grid_cfg:
                ib = grid_cfg["inner_box"]
                args.extend(
                    [
                        "--inner_box",
                        ", ".join(map(str, ib)) if isinstance(ib, list) else str(ib),
                    ]
                )
            if "outer_box" in grid_cfg:
                ob = grid_cfg["outer_box"]
                args.extend(
                    [
                        "--outer_box",
                        ", ".join(map(str, ob)) if isinstance(ob, list) else str(ob),
                    ]
                )
        lp_cfg = glide_cfg.get("ligprep", {})
        if lp_cfg:
            if "target_ph" in lp_cfg:
                args.extend(["--ligand_ph", str(lp_cfg["target_ph"])])
            if "ph_tolerance" in lp_cfg:
                args.extend(["--ligand_ph_tol", str(lp_cfg["ph_tolerance"])])
            if "force_field" in lp_cfg:
                args.extend(["--ligprep_ff", str(lp_cfg["force_field"])])
            if "max_states" in lp_cfg:
                args.extend(["--max_states", str(lp_cfg["max_states"])])

    # 6. Molecular Dynamics Simulation (Strict Validation)
    md_cfg = config.get("md_simulation", {})
    if md_cfg.get("enabled", False):
        args.extend(["--run_md", "true"])

        # System Build
        sys_build = md_cfg.get("system_build")
        if not sys_build:
            raise ValueError(
                "Config validation failed: 'md_simulation.system_build' section in md.yaml is required when MD is enabled!"
            )
        for key, flag in [
            ("water_model", "--water_model"),
            ("box_padding", "--box_padding"),
            ("salt_conc", "--md_salt_conc"),
        ]:
            if key not in sys_build:
                raise ValueError(
                    f"Config validation failed: 'md_simulation.system_build.{key}' in md.yaml is required!"
                )
            args.extend([flag, str(sys_build[key])])

        # Equilibration - STRICT VALIDATION (No hidden defaults)
        eq_cfg = md_cfg.get("equilibration")
        if not eq_cfg:
            raise ValueError(
                "Config validation failed: 'md_simulation.equilibration' section in md.yaml is required when MD is enabled!"
            )
        for key in ["nvt_short_ps", "nvt_long_ps", "npt_berendsen_ps", "npt_pr_ps"]:
            if key not in eq_cfg:
                raise ValueError(
                    f"Config validation failed: 'md_simulation.equilibration.{key}' in md.yaml is required!"
                )
            args.extend([f"--{key}", str(eq_cfg[key])])

        # Production - STRICT VALIDATION
        prod_cfg = md_cfg.get("production")
        if not prod_cfg:
            raise ValueError(
                "Config validation failed: 'md_simulation.production' section in md.yaml is required when MD is enabled!"
            )
        for key, flag in [
            ("total_time_ns", "--md_total_ns"),
            ("chunk_size_ns", "--md_chunk_ns"),
            ("replicates", "--md_replicates"),
            ("temperature", "--md_temperature"),
        ]:
            if key not in prod_cfg:
                raise ValueError(
                    f"Config validation failed: 'md_simulation.production.{key}' in md.yaml is required!"
                )
            args.extend([flag, str(prod_cfg[key])])

        # Selection & Trajectory
        if "pbc_center" in md_cfg:
            args.extend(["--pbc_center", str(md_cfg["pbc_center"])])
        if "pbc_extract" in md_cfg:
            args.extend(["--pbc_extract", str(md_cfg["pbc_extract"])])
        sel_cfg = md_cfg.get("selection", {})
        if sel_cfg.get("top_n"):
            args.extend(["--top_n_hits", str(sel_cfg["top_n"])])
        if sel_cfg.get("selected_hits"):
            hits = sel_cfg["selected_hits"]
            args.extend(
                ["--md_hits", ",".join(hits) if isinstance(hits, list) else str(hits)]
            )
    else:
        args.extend(["--run_md", "false"])

    # 7. MM-PBSA
    mmpbsa_cfg = config.get("mmpbsa", {})
    if mmpbsa_cfg.get("enabled", False):
        if "start_ns" in mmpbsa_cfg:
            args.extend(["--mmpbsa_start_ns", str(mmpbsa_cfg["start_ns"])])
        if "step_frames" in mmpbsa_cfg:
            args.extend(["--mmpbsa_step_frames", str(mmpbsa_cfg["step_frames"])])
        if "per_residue_decomp" in mmpbsa_cfg:
            args.extend(
                [
                    "--per_residue_decomp",
                    "true" if mmpbsa_cfg["per_residue_decomp"] else "false",
                ]
            )

    # Output directory
    args.extend(["--outdir", str(project_dir.resolve())])

    return args
