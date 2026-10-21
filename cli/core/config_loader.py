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
    expected_configs = [
        "project.yaml",
        "glide.yaml",
        "gnina.yaml",
        "md.yaml",
        "membrane.yaml",
        "mmpbsa.yaml",
        "manifest.yaml",
    ]

    for sdir in search_dirs:
        for cfg_name in expected_configs:
            yf = sdir / cfg_name
            if yf.exists():
                try:
                    with open(yf, "r", encoding="utf-8") as f:
                        data = yaml.safe_load(f)
                        if isinstance(data, dict):
                            merged.update(data)
                except Exception as e:
                    raise ValueError(
                        f"FATAL: Failed to parse required YAML file {yf.name}: {e}"
                    )

    return merged


def _build_prep_args(config: Dict[str, Any]) -> list:
    args = []
    prep_cfg = config.get("preparation")
    if not prep_cfg:
        raise ValueError("FATAL: 'preparation' config section is missing!")

    # Explicit mapping: (yaml_key, nextflow_flag)
    prep_mappings = [
        ("target_ph", "target_ph"),
        ("ph_tolerance", "ph_tolerance"),
        ("minimize_rmsd", "minimize_rmsd"),
    ]
    for yaml_k, nf_k in prep_mappings:
        if yaml_k not in prep_cfg:
            raise ValueError(f"FATAL: 'preparation.{yaml_k}' is required!")
        args.extend([f"--{nf_k}", str(prep_cfg[yaml_k])])

    if "cap_termini" in prep_cfg:
        args.extend(["--cap_termini", "true" if prep_cfg["cap_termini"] else "false"])

    # LigPrep parameters for chemical library and native ligand preparation
    glide_cfg = config.get("glide", {})
    lp_cfg = glide_cfg.get("ligprep", {})
    lig_ph = lp_cfg.get("target_ph", prep_cfg.get("target_ph", 7.4))
    lig_tol = lp_cfg.get("ph_tolerance", prep_cfg.get("ph_tolerance", 2.0))
    lig_ff = lp_cfg.get("force_field", 16)
    max_st = lp_cfg.get("max_states", 32)
    gnina_cfg = config.get("gnina", {})
    b_size = gnina_cfg.get("batch_size", 10)

    args.extend(
        [
            "--ligand_ph",
            str(lig_ph),
            "--ligand_ph_tol",
            str(lig_tol),
            "--ligprep_ff",
            str(lig_ff),
            "--max_states",
            str(max_st),
            "--batch_size",
            str(b_size),
        ]
    )

    mem_cfg = config.get("membrane", {})
    if mem_cfg.get("enabled"):
        args.extend(["--membrane", "true"])
        mem_mappings = [("type", "membrane_type"), ("topology", "membrane_topo")]
        for yaml_k, nf_k in mem_mappings:
            if yaml_k not in mem_cfg:
                raise ValueError(
                    f"FATAL: 'membrane.{yaml_k}' is required when membrane is enabled!"
                )
            args.extend([f"--{nf_k}", str(mem_cfg[yaml_k])])

        if "chains" in mem_cfg:
            chains = mem_cfg["chains"]
            args.extend(
                [
                    "--membrane_chains",
                    ",".join(chains) if isinstance(chains, list) else str(chains),
                ]
            )
    else:
        args.extend(["--membrane", "false"])

    return args


def _build_screen_args(config: Dict[str, Any], project_dir: Path) -> list:
    args = []
    screen_cfg = config.get("screening", {})
    if "library_file" in screen_cfg:
        lib_path = project_dir / screen_cfg["library_file"]
        if not lib_path.exists():
            lib_path = project_dir / "ligands" / screen_cfg["library_file"]
        if lib_path.exists():
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

    glide_cfg = config.get("glide", {})
    gnina_cfg = config.get("gnina", {})
    g_en = glide_cfg.get("enabled", False)
    n_en = gnina_cfg.get("enabled", False)

    if g_en and n_en:
        args.extend(["--dock_engine", "both"])
    elif g_en:
        args.extend(["--dock_engine", "glide"])
    elif n_en:
        args.extend(["--dock_engine", "gnina"])
    else:
        raise ValueError(
            "FATAL: At least one docking engine (glide or gnina) must be enabled for screening!"
        )

    if n_en:
        gnina_mappings = [
            ("exhaustiveness", "exhaustiveness"),
            ("cnn_scoring", "cnn_scoring"),
            ("num_modes", "num_modes"),
            ("autobox_add", "autobox_add"),
        ]
        for yaml_k, nf_k in gnina_mappings:
            if yaml_k not in gnina_cfg:
                raise ValueError(
                    f"FATAL: 'gnina.{yaml_k}' is required when GNINA is enabled!"
                )
            args.extend([f"--{nf_k}", str(gnina_cfg[yaml_k])])

        if "cnn_model" in gnina_cfg:
            args.extend(["--cnn_model", str(gnina_cfg["cnn_model"])])
        if "batch_size" in gnina_cfg:
            args.extend(["--batch_size", str(gnina_cfg["batch_size"])])

    if g_en:
        glide_mappings = [
            ("precision", "glide_precision"),
            ("poses_per_lig", "poses_per_lig"),
            ("forcefield", "glide_forcefield"),
            ("recep_vscale", "glide_recep_vscale"),
            ("lig_vscale", "glide_lig_vscale"),
            ("recep_ccut", "glide_recep_ccut"),
            ("lig_ccut", "glide_lig_ccut"),
        ]
        for yaml_k, nf_k in glide_mappings:
            if yaml_k not in glide_cfg:
                raise ValueError(
                    f"FATAL: 'glide.{yaml_k}' is required when Glide is enabled!"
                )
            args.extend([f"--{nf_k}", str(glide_cfg[yaml_k])])

        grid_cfg = glide_cfg.get("grid")
        if not grid_cfg:
            raise ValueError("FATAL: 'glide.grid' section is required!")

        for yaml_k in ["inner_box", "outer_box"]:
            if yaml_k not in grid_cfg:
                raise ValueError(f"FATAL: 'glide.grid.{yaml_k}' is required!")
            val = grid_cfg[yaml_k]
            args.extend(
                [
                    f"--{yaml_k}",
                    ", ".join(map(str, val)) if isinstance(val, list) else str(val),
                ]
            )

        lp_cfg = glide_cfg.get("ligprep", {})
        if lp_cfg:
            for yaml_k, nf_k in [
                ("target_ph", "ligand_ph"),
                ("ph_tolerance", "ligand_ph_tol"),
                ("force_field", "ligprep_ff"),
                ("max_states", "max_states"),
            ]:
                if yaml_k in lp_cfg:
                    args.extend([f"--{nf_k}", str(lp_cfg[yaml_k])])

    return args


def _build_rescore_args(config: Dict[str, Any]) -> list:
    args = []
    mem_cfg = config.get("membrane", {})
    if mem_cfg.get("enabled"):
        args.extend(["--membrane", "true"])
    else:
        args.extend(["--membrane", "false"])
    return args


def _build_md_setup_args(config: Dict[str, Any]) -> list:
    args = []
    md_cfg = config.get("md_simulation")
    if not md_cfg:
        raise ValueError("FATAL: 'md_simulation' section missing in md.yaml!")

    sys_build = md_cfg.get("system_build")
    if not sys_build:
        raise ValueError("FATAL: 'md_simulation.system_build' is required!")

    sys_mappings = [
        ("protein_forcefield", "protein_forcefield"),
        ("gaff_forcefield", "gaff_forcefield"),
        ("lipid_forcefield", "lipid_forcefield"),
        ("water_leaprc", "water_leaprc"),
        ("water_box", "water_box"),
        ("box_type", "box_type"),
        ("box_padding", "box_padding"),
        ("ss_cutoff", "ss_cutoff"),
        ("salt_conc", "md_salt_conc"),
    ]
    for yaml_k, nf_k in sys_mappings:
        if yaml_k not in sys_build:
            raise ValueError(
                f"FATAL: 'md_simulation.system_build.{yaml_k}' is required!"
            )
        args.extend([f"--{nf_k}", str(sys_build[yaml_k])])

    args.extend(["--water_names", ",".join(sys_build.get("water_names", []))])
    args.extend(["--ion_names", ",".join(sys_build.get("ion_names", []))])

    mem_cfg = config.get("membrane", {})
    if mem_cfg.get("enabled"):
        args.extend(["--membrane", "true"])
        lipids = mem_cfg.get("lipids", {})
        upper = lipids.get("upper_leaflet", {})
        lower = lipids.get("lower_leaflet", {})
        for k, v in [
            ("--upper_lipids", upper.get("composition")),
            ("--upper_ratio", upper.get("ratio")),
            ("--lower_lipids", lower.get("composition")),
            ("--lower_ratio", lower.get("ratio")),
        ]:
            if not v:
                raise ValueError(
                    f"FATAL: Membrane is enabled but {k.replace('--','')} is missing in membrane.yaml!"
                )
            args.extend([k, str(v)])

        pad_cfg = mem_cfg.get("padding", {})
        if "lipid_padding" in pad_cfg:
            args.extend(["--lipid_padding", str(pad_cfg["lipid_padding"])])
        args.extend(["--lipid_names", ",".join(mem_cfg.get("lipid_names", []))])
    else:
        args.extend(["--membrane", "false"])
        args.extend(["--lipid_names", ""])

    declash = md_cfg.get("declash")
    if not declash:
        raise ValueError("FATAL: 'md_simulation.declash' is required!")

    declash_mappings = [
        ("cutoff", "declash_cutoff"),
        ("niter", "declash_niter"),
        ("damping", "declash_damping"),
        ("max_move", "declash_max_move"),
        ("jitter", "declash_jitter"),
        ("tolerance", "declash_tol"),
    ]
    for yaml_k, nf_k in declash_mappings:
        if yaml_k not in declash:
            raise ValueError(f"FATAL: 'md_simulation.declash.{yaml_k}' is required!")
        args.extend([f"--{nf_k}", str(declash[yaml_k])])

    rf_list = sys_build.get("restraint_forces")
    if not isinstance(rf_list, list) or not rf_list:
        raise ValueError(
            "FATAL: 'restraint_forces' must be a non-empty list of force constants!"
        )
    args.extend(["--md_restraints_fc", ",".join(map(str, rf_list))])

    prod_cfg = md_cfg.get("production")
    if (
        prod_cfg
        and "restraints_fc" in prod_cfg
        and int(float(prod_cfg["restraints_fc"])) > 0
    ):
        args.extend(
            ["--prod_restraints_fc", str(int(float(prod_cfg["restraints_fc"])))]
        )
        args.extend(
            [
                "--prod_restraints_group",
                str(prod_cfg.get("restraints_group", "C-alpha")),
            ]
        )

    return args


def _build_dynamics_args(config: Dict[str, Any], project_dir: Path) -> list:
    args = ["--run_md", "true"]
    md_cfg = config.get("md_simulation")
    if not md_cfg:
        raise ValueError("FATAL: 'md_simulation' section missing in md.yaml!")

    exec_cfg = md_cfg.get("execution")
    if not exec_cfg:
        raise ValueError("FATAL: 'md_simulation.execution' block is required!")
    for k in ["minim_flags", "eq_flags", "prod_flags"]:
        if k not in exec_cfg:
            raise ValueError(f"FATAL: 'md_simulation.execution.{k}' is required!")
        args.extend([f"--{k}", str(exec_cfg[k])])

    prod_cfg = md_cfg.get("production")
    if not prod_cfg:
        raise ValueError("FATAL: 'md_simulation.production' is required!")

    prod_mappings = [
        ("total_time_ns", "md_total_ns"),
        ("chunk_size_ns", "md_chunk_ns"),
        ("replicates", "md_replicates"),
        ("dt_ps", "md_dt"),
        ("compressibility", "md_compressibility"),
        ("tau_p", "md_tau_p"),
    ]
    for yaml_k, nf_k in prod_mappings:
        if yaml_k not in prod_cfg:
            raise ValueError(f"FATAL: 'md_simulation.production.{yaml_k}' is required!")
        args.extend([f"--{nf_k}", str(prod_cfg[yaml_k])])

    if "restraints_fc" in prod_cfg and int(float(prod_cfg["restraints_fc"])) > 0:
        args.extend(
            ["--prod_restraints_fc", str(int(float(prod_cfg["restraints_fc"])))]
        )
        args.extend(
            [
                "--prod_restraints_group",
                str(prod_cfg.get("restraints_group", "C-alpha")),
            ]
        )

    protocol_dir = project_dir / "configs" / "md_protocol"
    if not protocol_dir.exists():
        raise FileNotFoundError(
            f"FATAL: MD Protocol directory not found: {protocol_dir}. Run 'sbdd md-init' first."
        )
    args.extend(["--md_protocol", str(protocol_dir.resolve())])

    return args


def _build_analyze_args(config: Dict[str, Any]) -> list:
    args = []
    md_cfg = config.get("md_simulation", {})
    metrics = md_cfg.get("metrics")
    if metrics:
        for k, flag in [
            ("heavy_cutoff", "--metric_heavy"),
            ("hyd_cutoff", "--metric_hyd"),
            ("pol_cutoff", "--metric_pol"),
        ]:
            if k in metrics:
                args.extend([flag, str(metrics[k])])

    if "pbc_center" in md_cfg:
        args.extend(["--pbc_center", str(md_cfg["pbc_center"])])
    if "pbc_extract" in md_cfg:
        args.extend(["--pbc_extract", str(md_cfg["pbc_extract"])])

    sys_build = md_cfg.get("system_build", {})
    args.extend(["--water_names", ",".join(sys_build.get("water_names", []))])
    args.extend(["--ion_names", ",".join(sys_build.get("ion_names", []))])

    mem_cfg = config.get("membrane", {})
    if mem_cfg.get("enabled"):
        args.extend(["--membrane", "true"])
        args.extend(["--lipid_names", ",".join(mem_cfg.get("lipid_names", []))])
        if "thickness" in mem_cfg:
            args.extend(["--membrane_thickness", str(mem_cfg["thickness"])])
    else:
        args.extend(["--membrane", "false"])
        args.extend(["--lipid_names", ""])

    mmpbsa_cfg = config.get("mmpbsa", {})
    if mmpbsa_cfg.get("enabled"):
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

    prod = md_cfg.get("production", {})
    temp = prod.get("temperature") or mmpbsa_cfg.get("temperature")
    salt = sys_build.get("salt_conc") or mmpbsa_cfg.get("salt_conc")
    if temp:
        args.extend(["--md_temperature", str(temp)])
    if salt:
        args.extend(["--md_salt_conc", str(salt)])

    return args


def validate_ligand_parametrization_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """Strict Fail-Fast validation for small molecule charge and forcefield parametrization."""
    md_cfg = config.get("md_simulation")
    if not md_cfg:
        raise ValueError("FATAL: 'md_simulation' section missing in md.yaml!")

    param_cfg = md_cfg.get("ligand_parametrization")
    if not param_cfg:
        raise ValueError(
            "FATAL: 'md_simulation.ligand_parametrization' section missing in md.yaml! "
            "Parametrization method ('resp' or 'bcc') and target forcefield must be explicitly declared."
        )

    method = param_cfg.get("method")
    if not method or method not in ["resp", "bcc"]:
        raise ValueError(
            f"FATAL: Invalid or missing 'md_simulation.ligand_parametrization.method': '{method}'. "
            "Supported methods are 'resp' (ORCA 6 quantum ESP) or 'bcc' (AM1-BCC)."
        )

    ff = param_cfg.get("forcefield")
    if not ff or ff not in ["gaff", "gaff2"]:
        raise ValueError(
            f"FATAL: Invalid or missing 'md_simulation.ligand_parametrization.forcefield': '{ff}'. "
            "Supported small molecule force fields are 'gaff' or 'gaff2'."
        )

    if method == "resp":
        orca_cfg = param_cfg.get("orca")
        if not orca_cfg or not isinstance(orca_cfg, dict):
            raise ValueError(
                "FATAL: 'md_simulation.ligand_parametrization.orca' section is required when method='resp'!"
            )
        required_orca_keys = [
            "nprocs",
            "maxcore_mb",
            "opt_theory",
            "opt_basis",
            "opt_dispersion",
            "opt_solvent",
            "opt_rijcosx",
            "opt_maxiter",
            "scf_convergence",
            "scf_maxiter",
            "esp_theory",
            "esp_basis",
            "esp_density",
        ]
        for k in required_orca_keys:
            if k not in orca_cfg or orca_cfg[k] is None or orca_cfg[k] == "":
                raise ValueError(
                    f"FATAL: Required key 'md_simulation.ligand_parametrization.orca.{k}' is missing or empty in md.yaml!"
                )
        if int(orca_cfg["nprocs"]) <= 0:
            raise ValueError(
                f"FATAL: orca.nprocs must be a positive integer, got: {orca_cfg['nprocs']}"
            )
        if int(orca_cfg["maxcore_mb"]) <= 0:
            raise ValueError(
                f"FATAL: orca.maxcore_mb must be a positive integer, got: {orca_cfg['maxcore_mb']}"
            )
        if int(orca_cfg["opt_maxiter"]) <= 0:
            raise ValueError(
                f"FATAL: orca.opt_maxiter must be a positive integer, got: {orca_cfg['opt_maxiter']}"
            )
        if int(orca_cfg["scf_maxiter"]) <= 0:
            raise ValueError(
                f"FATAL: orca.scf_maxiter must be a positive integer, got: {orca_cfg['scf_maxiter']}"
            )

    elif method == "bcc":
        bcc_cfg = param_cfg.get("bcc")
        if not bcc_cfg or not isinstance(bcc_cfg, dict):
            raise ValueError(
                "FATAL: 'md_simulation.ligand_parametrization.bcc' section is required when method='bcc'!"
            )
        required_bcc_keys = ["scfconv", "grms_tol"]
        for k in required_bcc_keys:
            if k not in bcc_cfg or bcc_cfg[k] is None or bcc_cfg[k] == "":
                raise ValueError(
                    f"FATAL: Required key 'md_simulation.ligand_parametrization.bcc.{k}' is missing or empty in md.yaml!"
                )

    return param_cfg


def _build_parametrize_args(config: Dict[str, Any]) -> list:
    """Builds Nextflow command line arguments for the parametrization workflow."""
    param_cfg = validate_ligand_parametrization_config(config)
    method = param_cfg["method"]
    ff = param_cfg["forcefield"]
    args = [
        "--param_method",
        method,
        "--param_ff",
        ff,
    ]
    if method == "resp":
        orca = param_cfg["orca"]
        args.extend(
            [
                "--orca_nprocs",
                str(orca["nprocs"]),
                "--orca_maxcore",
                str(orca["maxcore_mb"]),
                "--orca_opt_theory",
                str(orca["opt_theory"]),
                "--orca_opt_basis",
                str(orca["opt_basis"]),
                "--orca_opt_dispersion",
                str(orca["opt_dispersion"]),
                "--orca_opt_solvent",
                str(orca["opt_solvent"]),
                "--orca_opt_rijcosx",
                "true" if orca["opt_rijcosx"] else "false",
                "--orca_opt_maxiter",
                str(orca["opt_maxiter"]),
                "--orca_scf_conv",
                str(orca["scf_convergence"]),
                "--orca_scf_maxiter",
                str(orca["scf_maxiter"]),
                "--orca_esp_theory",
                str(orca["esp_theory"]),
                "--orca_esp_basis",
                str(orca["esp_basis"]),
                "--orca_esp_density",
                str(orca["esp_density"]),
            ]
        )
    elif method == "bcc":
        bcc = param_cfg["bcc"]
        args.extend(
            [
                "--bcc_scfconv",
                str(bcc["scfconv"]),
                "--bcc_grms_tol",
                str(bcc["grms_tol"]),
            ]
        )
    return args


def build_nextflow_args(
    config: Dict[str, Any], project_dir: Path, entrypoint: str
) -> list:
    """
    Translates dictionary of parameters into Nextflow CLI arguments with strict routing validation.
    Zero cross-contamination between domains.
    """
    args = []

    # 1. Base Project Validations (Universal)
    proj_cfg = config.get("project", {})
    if not proj_cfg.get("name"):
        raise ValueError("FATAL: 'project.name' in project.yaml is required!")

    args.extend(
        [
            "--project_name",
            str(proj_cfg["name"]),
            "--outdir",
            str(project_dir.resolve()),
        ]
    )

    # 2. Strict Domain Routing
    if entrypoint == "prepare":
        args.extend(_build_prep_args(config))
    elif entrypoint in ["screen", "dock"]:
        args.extend(_build_screen_args(config, project_dir))
    elif entrypoint == "rescore":
        args.extend(_build_rescore_args(config))
    elif entrypoint == "md-setup":
        args.extend(_build_md_setup_args(config))
    elif entrypoint == "dynamics":
        args.extend(_build_dynamics_args(config, project_dir))
    elif entrypoint == "md-analyze":
        args.extend(_build_analyze_args(config))
    elif entrypoint == "heal":
        pass  # Heal only requires the base project params mapped above
    elif entrypoint == "parametrize":
        args.extend(_build_parametrize_args(config))
    else:
        raise ValueError(
            f"FATAL: Unknown routing entrypoint '{entrypoint}'. Cannot build arguments."
        )

    return args
