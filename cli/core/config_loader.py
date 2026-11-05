"""
File: cli/core/config_loader.py
Description: Scans, loads, and strictly validates atomic configuration YAML files from configs/.
Zero silent fallbacks: required parameters must be explicitly specified.
"""

from pathlib import Path
import re
from typing import Dict, Any
import yaml

SLURM_GRES_CONF_PATH = Path("/etc/slurm/gres.conf")


def _resolve_cluster_gpu_shards(
    gres_conf_path: Path = SLURM_GRES_CONF_PATH,
) -> int:
    """
    Reads total GPU shard count per device directly from the cluster's /etc/slurm/gres.conf.
    Enforces strict Fail-Fast validation with zero hardcoded fallbacks.
    """
    if not gres_conf_path.exists():
        raise FileNotFoundError(
            f"FATAL: Slurm GRES configuration file not found at {gres_conf_path}! "
            "Cannot determine cluster GPU shard capacity."
        )

    content = gres_conf_path.read_text(encoding="utf-8")
    for raw_line in content.splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if re.search(r"\bName\s*=\s*shard\b", line, re.IGNORECASE):
            m_count = re.search(r"\bCount\s*=\s*(\d+)\b", line, re.IGNORECASE)
            if not m_count:
                raise ValueError(
                    f"FATAL: Malformed 'Name=shard' entry without 'Count=<N>' in {gres_conf_path}: '{raw_line}'"
                )
            total_shards = int(m_count.group(1))
            if total_shards <= 0:
                raise ValueError(
                    f"FATAL: Invalid non-positive shard Count={total_shards} in {gres_conf_path}!"
                )
            return total_shards

    raise ValueError(
        f"FATAL: No 'Name=shard Count=<N>' declaration found in {gres_conf_path}!"
    )


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

    # LigPrep parameters for chemical library and native ligand preparation (Strict Fail-Fast)
    glide_cfg = config.get("glide")
    if not glide_cfg or "ligprep" not in glide_cfg:
        raise ValueError("FATAL: 'glide.ligprep' section is required in glide.yaml!")
    lp_cfg = glide_cfg["ligprep"]
    for req_k in [
        "target_ph",
        "ph_tolerance",
        "force_field",
        "max_states",
        "batch_size",
        "cpus",
        "memory_gb",
    ]:
        if req_k not in lp_cfg or lp_cfg[req_k] is None:
            raise ValueError(
                f"FATAL: 'glide.ligprep.{req_k}' is required in glide.yaml!"
            )
    if (
        int(lp_cfg["batch_size"]) <= 0
        or int(lp_cfg["cpus"]) <= 0
        or int(lp_cfg["memory_gb"]) <= 0
    ):
        raise ValueError(
            "FATAL: 'glide.ligprep.batch_size', 'cpus', and 'memory_gb' must be positive integers!"
        )

    lig_ph = lp_cfg["target_ph"]
    lig_tol = lp_cfg["ph_tolerance"]
    lig_ff = lp_cfg["force_field"]
    max_st = lp_cfg["max_states"]
    b_size = int(lp_cfg["batch_size"])

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
            "--ligprep_cpus",
            str(int(lp_cfg["cpus"])),
            "--ligprep_mem_gb",
            str(int(lp_cfg["memory_gb"])),
        ]
    )

    mem_cfg = config.get("membrane", {})
    if mem_cfg.get("enabled"):
        args.extend(["--membrane", "true"])

        valid_ppm_codes = {
            "PMm", "PMp", "PMf", "Erf", "ERm", "GOL", "LYS", "END", "VAC",
            "MOM", "MIM", "THp", "THb", "GnO", "GnI", "GpI", "ARC",
            "LPC", "MPC", "OPC", "EPC", "MIC"
        }
        mem_type = str(mem_cfg.get("type", "")).strip()
        if mem_type not in valid_ppm_codes:
            raise ValueError(
                f"FATAL CONFIGURATION ERROR: Invalid 'membrane.type': '{mem_type}'! "
                f"PPM 3.0 requires an exact 3-letter code from the official table (case-sensitive):\n"
                f"{sorted(list(valid_ppm_codes))}"
            )

        mem_topo = str(mem_cfg.get("topology", "")).strip().lower()
        if mem_topo not in ["in", "out"]:
            raise ValueError(
                f"FATAL CONFIGURATION ERROR: Invalid 'membrane.topology': '{mem_topo}'! Must be 'in' or 'out'."
            )

        args.extend(["--membrane_type", mem_type, "--membrane_topo", mem_topo])

        if "chains" not in mem_cfg or not mem_cfg["chains"]:
            raise ValueError(
                "FATAL: 'membrane.chains' is required in membrane.yaml when membrane is enabled!"
            )
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
            ("batch_size", "gnina_batch_size"),
        ]
        for yaml_k, nf_k in gnina_mappings:
            if yaml_k not in gnina_cfg or gnina_cfg[yaml_k] is None:
                raise ValueError(
                    f"FATAL: 'gnina.{yaml_k}' is required in gnina.yaml when GNINA is enabled!"
                )
            args.extend([f"--{nf_k}", str(gnina_cfg[yaml_k])])
        if int(gnina_cfg["batch_size"]) <= 0:
            raise ValueError("FATAL: 'gnina.batch_size' must be a positive integer!")

        g_res = gnina_cfg.get("resources")
        if not isinstance(g_res, dict):
            raise ValueError(
                "FATAL: 'gnina.resources' block is required in gnina.yaml!"
            )
        for rk in ["cpus", "memory_gb", "slurm_gres"]:
            if rk not in g_res or g_res[rk] is None or str(g_res[rk]).strip() == "":
                raise ValueError(
                    f"FATAL: 'gnina.resources.{rk}' is required in gnina.yaml!"
                )
        g_cpus = int(g_res["cpus"])
        g_mem = int(g_res["memory_gb"])
        g_gres = str(g_res["slurm_gres"]).strip()
        if g_cpus <= 0 or g_mem <= 0:
            raise ValueError(
                "FATAL: 'gnina.resources.cpus' and 'memory_gb' must be positive integers!"
            )
        args.extend(
            [
                "--gnina_cpus",
                str(g_cpus),
                "--gnina_mem_gb",
                str(g_mem),
                "--gnina_gres",
                g_gres,
            ]
        )

        if (
            "cnn_model" in gnina_cfg
            and gnina_cfg["cnn_model"] is not None
            and str(gnina_cfg["cnn_model"]).strip() != ""
        ):
            args.extend(["--cnn_model", str(gnina_cfg["cnn_model"]).strip()])
    else:
        args.extend(["--gnina_batch_size", "0"])

    if g_en:
        glide_mappings = [
            ("precision", "glide_precision"),
            ("poses_per_lig", "poses_per_lig"),
            ("forcefield", "glide_forcefield"),
            ("recep_vscale", "glide_recep_vscale"),
            ("lig_vscale", "glide_lig_vscale"),
            ("recep_ccut", "glide_recep_ccut"),
            ("lig_ccut", "glide_lig_ccut"),
            ("batch_size", "glide_batch_size"),
            ("postdock", "glide_postdock"),
        ]
        for yaml_k, nf_k in glide_mappings:
            if yaml_k not in glide_cfg or glide_cfg[yaml_k] is None:
                raise ValueError(
                    f"FATAL: 'glide.{yaml_k}' is required in glide.yaml when Glide is enabled!"
                )
            args.extend([f"--{nf_k}", str(glide_cfg[yaml_k])])
        if int(glide_cfg["batch_size"]) <= 0:
            raise ValueError("FATAL: 'glide.batch_size' must be a positive integer!")

        gl_res = glide_cfg.get("resources")
        if not isinstance(gl_res, dict):
            raise ValueError(
                "FATAL: 'glide.resources' block is required in glide.yaml!"
            )
        for rk in ["cpus", "memory_gb"]:
            if rk not in gl_res or gl_res[rk] is None:
                raise ValueError(
                    f"FATAL: 'glide.resources.{rk}' is required in glide.yaml!"
                )
        if int(gl_res["cpus"]) <= 0 or int(gl_res["memory_gb"]) <= 0:
            raise ValueError(
                "FATAL: 'glide.resources.cpus' and 'memory_gb' must be positive integers!"
            )
        args.extend(
            [
                "--glide_cpus",
                str(int(gl_res["cpus"])),
                "--glide_mem_gb",
                str(int(gl_res["memory_gb"])),
            ]
        )

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
    else:
        args.extend(["--glide_batch_size", "0"])

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
        ("hmr_mass", "hmr_mass"),
    ]
    for yaml_k, nf_k in sys_mappings:
        if yaml_k not in sys_build or sys_build[yaml_k] is None:
            raise ValueError(
                f"FATAL: 'md_simulation.system_build.{yaml_k}' is required!"
            )
        args.extend([f"--{nf_k}", str(sys_build[yaml_k])])

    if "hmr" not in sys_build or sys_build["hmr"] is None:
        raise ValueError("FATAL: 'md_simulation.system_build.hmr' boolean is required in md.yaml!")
    if float(sys_build["hmr_mass"]) <= 1.008:
        raise ValueError("FATAL: 'md_simulation.system_build.hmr_mass' must be > 1.008 Da!")
    args.extend(["--hmr", "true" if sys_build["hmr"] else "false"])

    water_names = sys_build.get("water_names")
    if not isinstance(water_names, list) or not water_names:
        raise ValueError(
            "FATAL: 'md_simulation.system_build.water_names' must be a non-empty list in md.yaml!"
        )
    ion_names = sys_build.get("ion_names")
    if not isinstance(ion_names, list) or not ion_names:
        raise ValueError(
            "FATAL: 'md_simulation.system_build.ion_names' must be a non-empty list in md.yaml!"
        )
    args.extend(["--water_names", ",".join(water_names)])
    args.extend(["--ion_names", ",".join(ion_names)])

    mem_cfg = config.get("membrane")
    if not isinstance(mem_cfg, dict) or "enabled" not in mem_cfg:
        raise ValueError(
            "FATAL: 'membrane' section with explicit 'enabled' boolean is required in membrane.yaml!"
        )

    if mem_cfg["enabled"]:
        args.extend(["--membrane", "true"])
        lipids = mem_cfg.get("lipids")
        if not isinstance(lipids, dict):
            raise ValueError(
                "FATAL: 'membrane.lipids' section is required in membrane.yaml when membrane is enabled!"
            )
        upper = lipids.get("upper_leaflet")
        lower = lipids.get("lower_leaflet")
        if not isinstance(upper, dict) or not isinstance(lower, dict):
            raise ValueError(
                "FATAL: 'membrane.lipids.upper_leaflet' and 'lower_leaflet' are required in membrane.yaml!"
            )
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

        pad_cfg = mem_cfg.get("padding")
        if (
            not isinstance(pad_cfg, dict)
            or "lipid_padding" not in pad_cfg
            or pad_cfg["lipid_padding"] is None
            or "box_padding" not in pad_cfg
            or pad_cfg["box_padding"] is None
        ):
            raise ValueError(
                "FATAL: 'membrane.padding.lipid_padding' and 'membrane.padding.box_padding' are required in membrane.yaml!"
            )
        args.extend(["--lipid_padding", str(pad_cfg["lipid_padding"])])
        # Override box_padding with membrane water layer thickness when membrane is enabled
        idx_bp = args.index("--box_padding")
        args[idx_bp + 1] = str(pad_cfg["box_padding"])

        pack_cfg = mem_cfg.get("packing")
        if (
            not isinstance(pack_cfg, dict)
            or "nloop" not in pack_cfg
            or pack_cfg["nloop"] is None
            or "nloop_all" not in pack_cfg
            or pack_cfg["nloop_all"] is None
        ):
            raise ValueError(
                "FATAL: 'membrane.packing.nloop' and 'membrane.packing.nloop_all' must be explicitly specified in membrane.yaml!"
            )
        if int(pack_cfg["nloop"]) <= 0 or int(pack_cfg["nloop_all"]) <= 0:
            raise ValueError(
                "FATAL: 'membrane.packing.nloop' and 'nloop_all' must be positive integers!"
            )
        args.extend(
            [
                "--packmol_nloop",
                str(int(pack_cfg["nloop"])),
                "--packmol_nloop_all",
                str(int(pack_cfg["nloop_all"])),
            ]
        )

        lipid_names = mem_cfg.get("lipid_names")
        if not isinstance(lipid_names, list) or not lipid_names:
            raise ValueError(
                "FATAL: 'membrane.lipid_names' must be a non-empty list in membrane.yaml!"
            )
        args.extend(["--lipid_names", ",".join(lipid_names)])
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
        if yaml_k not in declash or declash[yaml_k] is None:
            raise ValueError(f"FATAL: 'md_simulation.declash.{yaml_k}' is required!")
        args.extend([f"--{nf_k}", str(declash[yaml_k])])

    rf_list = sys_build.get("restraint_forces")
    if not isinstance(rf_list, list) or not rf_list:
        raise ValueError(
            "FATAL: 'restraint_forces' must be a non-empty list of force constants!"
        )
    args.extend(["--md_restraints_fc", ",".join(map(str, rf_list))])

    prod_cfg = md_cfg.get("production")
    if not isinstance(prod_cfg, dict) or "restraints_fc" not in prod_cfg or prod_cfg["restraints_fc"] is None:
        raise ValueError(
            "FATAL: 'md_simulation.production.restraints_fc' is required in md.yaml!"
        )
    prod_fc_val = int(float(prod_cfg["restraints_fc"]))
    args.extend(["--prod_restraints_fc", str(prod_fc_val)])
    if prod_fc_val > 0:
        if "restraints_group" not in prod_cfg or not prod_cfg["restraints_group"]:
            raise ValueError(
                "FATAL: 'md_simulation.production.restraints_group' is required when restraints_fc > 0!"
            )
        args.extend(["--prod_restraints_group", str(prod_cfg["restraints_group"])])

    return args


def _build_dynamics_args(config: Dict[str, Any], project_dir: Path) -> list:
    args = ["--run_md", "true"]
    md_cfg = config.get("md_simulation")
    if not md_cfg:
        raise ValueError("FATAL: 'md_simulation' section missing in md.yaml!")

    exec_cfg = md_cfg.get("execution")
    if not exec_cfg:
        raise ValueError("FATAL: 'md_simulation.execution' block is required!")
    for k in [
        "minim_flags",
        "eq_flags",
        "prod_flags",
        "cpus",
        "memory_gb",
        "jobs_per_gpu",
    ]:
        if k not in exec_cfg or exec_cfg[k] is None:
            raise ValueError(
                f"FATAL: 'md_simulation.execution.{k}' is required in md.yaml!"
            )
    md_cpus = int(exec_cfg["cpus"])
    md_mem_gb = int(exec_cfg["memory_gb"])
    jobs_per_gpu = int(exec_cfg["jobs_per_gpu"])
    if md_cpus <= 0 or md_mem_gb <= 0 or jobs_per_gpu <= 0:
        raise ValueError(
            "FATAL: 'md_simulation.execution.cpus', 'memory_gb', and 'jobs_per_gpu' must be positive integers!"
        )

    total_gpu_shards = _resolve_cluster_gpu_shards()
    if total_gpu_shards % jobs_per_gpu != 0:
        raise ValueError(
            f"FATAL: 'md_simulation.execution.jobs_per_gpu' ({jobs_per_gpu}) must evenly divide "
            f"the cluster GPU shard capacity ({total_gpu_shards} in {SLURM_GRES_CONF_PATH})!"
        )
    shards_per_job = total_gpu_shards // jobs_per_gpu

    for k in ["minim_flags", "eq_flags", "prod_flags"]:
        args.extend([f"--{k}", str(exec_cfg[k])])
    args.extend(
        [
            "--md_cpus",
            str(md_cpus),
            "--md_mem_gb",
            str(md_mem_gb),
            "--md_gres",
            f"shard:{shards_per_job}",
        ]
    )

    prod_cfg = md_cfg.get("production")
    if not prod_cfg:
        raise ValueError("FATAL: 'md_simulation.production' is required!")

    sys_build = md_cfg.get("system_build", {})
    if "hmr" not in sys_build or sys_build["hmr"] is None:
        raise ValueError("FATAL: 'md_simulation.system_build.hmr' is required in md.yaml!")

    prod_mappings = [
        ("total_time_ns", "md_total_ns"),
        ("chunk_size_ns", "md_chunk_ns"),
        ("replicates", "md_replicates"),
        ("dt_ps", "md_dt"),
        ("compressibility", "md_compressibility"),
        ("tau_p", "md_tau_p"),
    ]
    for yaml_k, nf_k in prod_mappings:
        if yaml_k not in prod_cfg or prod_cfg[yaml_k] is None:
            raise ValueError(f"FATAL: 'md_simulation.production.{yaml_k}' is required!")
        args.extend([f"--{nf_k}", str(prod_cfg[yaml_k])])

    if float(prod_cfg["dt_ps"]) > 0.0025 and not bool(sys_build["hmr"]):
        raise ValueError(
            "FATAL: Timestep dt_ps > 0.0025 ps requires 'md_simulation.system_build.hmr: true' in md.yaml!"
        )

    tot_ns = float(prod_cfg["total_time_ns"])
    chk_ns = float(prod_cfg["chunk_size_ns"])
    if chk_ns <= 0 or tot_ns <= 0:
        raise ValueError("FATAL: total_time_ns and chunk_size_ns must be positive!")
    num_chunks = int(round(tot_ns / chk_ns))
    if abs(tot_ns - num_chunks * chk_ns) > 1e-4:
        raise ValueError(
            f"FATAL: total_time_ns ({tot_ns}) must be an exact multiple of chunk_size_ns ({chk_ns})!"
        )
    if num_chunks < 1 or num_chunks > 20:
        raise ValueError(
            f"FATAL: Required number of MD chunks ({num_chunks}) exceeds supported static DAG limit (1..20)!"
        )

    if "restraints_fc" not in prod_cfg or prod_cfg["restraints_fc"] is None:
        raise ValueError(
            "FATAL: 'md_simulation.production.restraints_fc' is required in md.yaml!"
        )
    prod_fc_val = int(float(prod_cfg["restraints_fc"]))
    args.extend(["--prod_restraints_fc", str(prod_fc_val)])
    if prod_fc_val > 0:
        if "restraints_group" not in prod_cfg or not prod_cfg["restraints_group"]:
            raise ValueError(
                "FATAL: 'md_simulation.production.restraints_group' is required when restraints_fc > 0!"
            )
        args.extend(["--prod_restraints_group", str(prod_cfg["restraints_group"])])

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

    val_sampling = validate_analysis_sampling_config(config)
    args.extend(
        [
            "--analysis_start_ns",
            str(val_sampling["start_ns"]),
            "--analysis_sample_interval_ps",
            str(val_sampling["sample_interval_ps"]),
        ]
    )

    mmpbsa_cfg = config.get("mmpbsa", {})
    if mmpbsa_cfg.get("enabled"):
        val_mmpbsa = validate_mmpbsa_config(config)
        args.extend(
            [
                "--per_residue_decomp",
                "true" if val_mmpbsa["per_residue_decomp"] else "false",
                "--mmpbsa_method",
                str(val_mmpbsa["method"]),
                "--mmpbsa_gb_model",
                str(val_mmpbsa["gb_model"]),
                "--mmpbsa_emem",
                str(val_mmpbsa["emem"]),
                "--mmpbsa_indi",
                str(val_mmpbsa["indi"]),
                "--mmpbsa_fillratio",
                str(val_mmpbsa["pb_fillratio"]),
                "--mmpbsa_cpus",
                str(val_mmpbsa["cpus"]),
                "--mmpbsa_mem_gb",
                str(val_mmpbsa["memory_gb"]),
            ]
        )

    prod = md_cfg.get("production", {})
    temp = prod.get("temperature") or mmpbsa_cfg.get("temperature")
    salt = sys_build.get("salt_conc") or mmpbsa_cfg.get("salt_conc")
    if temp:
        args.extend(["--md_temperature", str(temp)])
    if salt:
        args.extend(["--md_salt_conc", str(salt)])

    # Validate and append clustering parameters if declared in configuration
    cluster_cfg = md_cfg.get("clustering")
    if cluster_cfg:
        val_cluster = validate_clustering_config(config)
        args.extend(
            [
                "--cluster_align_group",
                str(val_cluster["align_group"]),
                "--cluster_feature_group",
                str(val_cluster["cluster_group"]),
                "--cluster_dim_reduction",
                str(val_cluster["dim_reduction"]),
                "--cluster_algorithm",
                str(val_cluster["algorithm"]),
                "--cluster_rmsd_cutoff",
                str(val_cluster["rmsd_cutoff"]),
                "--cluster_tica_lag_ns",
                str(val_cluster["tica_lag_ns"]),
                "--cluster_variance_cutoff",
                str(val_cluster["variance_cutoff"]),
                "--cluster_n_clusters",
                str(val_cluster["n_clusters"]),
                "--cluster_pocket_radius",
                str(val_cluster["pocket_radius"]),
                "--cluster_min_pop",
                str(val_cluster["min_cluster_population"]),
                "--cluster_temp_k",
                str(val_cluster["temp_k"]),
                "--cluster_auto_promote",
                "true" if val_cluster["auto_promote"] else "false",
            ]
        )

    return args


def validate_analysis_sampling_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Central Single Source of Truth validation for trajectory sampling stride.
    Guarantees that MM-PBSA, Clustering, and Stratification share an identical time grid.
    """
    sampling_cfg = config.get("analysis_sampling") or config.get(
        "md_simulation", {}
    ).get("analysis_sampling")

    if not sampling_cfg or not isinstance(sampling_cfg, dict):
        # Legacy tolerance check: inspect mmpbsa section if present
        m_cfg = config.get("mmpbsa", {})
        if "start_ns" in m_cfg and "sample_interval_ps" in m_cfg:
            sampling_cfg = {
                "start_ns": m_cfg["start_ns"],
                "sample_interval_ps": m_cfg["sample_interval_ps"],
            }
        else:
            raise ValueError(
                "FATAL SSOT VIOLATION: Missing required 'analysis_sampling' configuration block in md.yaml!\n"
                "Both MM-PBSA and Clustering must derive their sampling stride from this single source of truth."
            )

    required_keys = ["start_ns", "sample_interval_ps"]
    for key in required_keys:
        if (
            key not in sampling_cfg
            or sampling_cfg[key] is None
            or sampling_cfg[key] == ""
        ):
            raise ValueError(
                f"FATAL: Missing required parameter 'analysis_sampling.{key}' in md.yaml! "
                "Silent fallbacks are strictly forbidden by architectural law."
            )

    start_ns = float(sampling_cfg["start_ns"])
    interval_ps = float(sampling_cfg["sample_interval_ps"])

    if start_ns < 0.0:
        raise ValueError(
            f"FATAL: 'analysis_sampling.start_ns' must be non-negative, got: {start_ns}"
        )
    if interval_ps <= 0.0:
        raise ValueError(
            f"FATAL: 'analysis_sampling.sample_interval_ps' must be positive, got: {interval_ps}"
        )

    return {"start_ns": start_ns, "sample_interval_ps": interval_ps}


def validate_mmpbsa_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """Strict Fail-Fast validation of MM-PBSA thermodynamic parameters."""
    mmpbsa_cfg = config.get("mmpbsa")
    if not mmpbsa_cfg or not isinstance(mmpbsa_cfg, dict):
        raise ValueError(
            "FATAL: 'mmpbsa' configuration block is missing or invalid in mmpbsa.yaml!"
        )

    # Inherit single source of truth for trajectory sampling
    sampling = validate_analysis_sampling_config(config)
    mmpbsa_cfg["start_ns"] = sampling["start_ns"]
    mmpbsa_cfg["sample_interval_ps"] = sampling["sample_interval_ps"]

    required_keys = ["temperature", "salt_conc", "per_residue_decomp"]
    for key in required_keys:
        if key not in mmpbsa_cfg or mmpbsa_cfg[key] is None or mmpbsa_cfg[key] == "":
            raise ValueError(
                f"FATAL: Missing required parameter 'mmpbsa.{key}' in mmpbsa.yaml! "
                "Silent fallbacks and default constants are strictly forbidden by architectural law."
            )

    if float(mmpbsa_cfg["temperature"]) <= 0.0:
        raise ValueError(
            f"FATAL: 'temperature' must be positive Kelvin, got: {mmpbsa_cfg['temperature']}"
        )

    if "method" not in mmpbsa_cfg or not mmpbsa_cfg["method"]:
        raise ValueError("FATAL: 'mmpbsa.method' is required in mmpbsa.yaml!")
    method = str(mmpbsa_cfg["method"]).strip().lower()
    if method not in ["auto", "gb", "pb"]:
        raise ValueError(
            f"FATAL: Invalid mmpbsa.method: '{method}'. Must be 'auto', 'gb', or 'pb'."
        )

    res_block = mmpbsa_cfg.get("resources")
    if not isinstance(res_block, dict):
        raise ValueError("FATAL: 'mmpbsa.resources' block is required in mmpbsa.yaml!")
    for r_key in ["cpus", "memory_gb"]:
        if r_key not in res_block or res_block[r_key] is None:
            raise ValueError(
                f"FATAL: 'mmpbsa.resources.{r_key}' is required in mmpbsa.yaml!"
            )
    if int(res_block["cpus"]) <= 0 or int(res_block["memory_gb"]) <= 0:
        raise ValueError(
            "FATAL: 'mmpbsa.resources.cpus' and 'memory_gb' must be positive integers!"
        )
    mmpbsa_cfg["cpus"] = int(res_block["cpus"])
    mmpbsa_cfg["memory_gb"] = int(res_block["memory_gb"])

    solv = mmpbsa_cfg.get("solvation")
    if not solv or not isinstance(solv, dict):
        raise ValueError("FATAL: 'mmpbsa.solvation' block is required in mmpbsa.yaml!")
    for s_key in ["gb_model", "pb_fillratio", "emem", "indi"]:
        if s_key not in solv or solv[s_key] is None:
            raise ValueError(
                f"FATAL: 'mmpbsa.solvation.{s_key}' is required in mmpbsa.yaml!"
            )

    mmpbsa_cfg["gb_model"] = int(solv["gb_model"])
    mmpbsa_cfg["pb_fillratio"] = float(solv["pb_fillratio"])
    mmpbsa_cfg["emem"] = float(solv["emem"])
    mmpbsa_cfg["indi"] = float(solv["indi"])

    return mmpbsa_cfg


def validate_clustering_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """Strict Fail-Fast validation of conformational landscape clustering parameters."""
    md_cfg = config.get("md_simulation")
    if not md_cfg:
        raise ValueError("FATAL: 'md_simulation' section missing in md.yaml!")

    cluster_cfg = md_cfg.get("clustering")
    if not cluster_cfg or not isinstance(cluster_cfg, dict):
        raise ValueError(
            "FATAL: 'md_simulation.clustering' configuration block is missing or invalid in md.yaml! "
            "All parameters must be explicitly specified."
        )

    # Inherit single source of truth for trajectory sampling
    sampling = validate_analysis_sampling_config(config)
    cluster_cfg["start_ns"] = sampling["start_ns"]
    cluster_cfg["sample_interval_ps"] = sampling["sample_interval_ps"]

    required_keys = [
        "align_group",
        "cluster_group",
        "dim_reduction",
        "algorithm",
        "tica_lag_ns",
        "variance_cutoff",
        "pocket_radius",
        "min_cluster_population",
        "temp_k",
        "auto_promote",
    ]

    for key in required_keys:
        if key not in cluster_cfg or cluster_cfg[key] is None or cluster_cfg[key] == "":
            raise ValueError(
                f"FATAL: Missing required clustering parameter 'md_simulation.clustering.{key}' in md.yaml! "
                "Silent fallbacks and default constants are strictly forbidden by architectural law."
            )

    dim_red = str(cluster_cfg["dim_reduction"]).strip().lower()
    if dim_red not in ["pca", "tica"]:
        raise ValueError(
            f"FATAL: Unsupported dimensionality reduction method '{dim_red}' in md.yaml! "
            "Must be strictly 'pca' or 'tica'."
        )

    algo = str(cluster_cfg["algorithm"]).strip().lower()
    if algo not in ["hdbscan", "hierarchical", "auto_kmeans", "kmeans"]:
        raise ValueError(
            f"FATAL: Unsupported clustering algorithm '{algo}' in md.yaml! "
            "Must be strictly 'hdbscan', 'hierarchical', 'auto_kmeans', or 'kmeans'."
        )

    for k_req in ["rmsd_cutoff", "n_clusters"]:
        if k_req not in cluster_cfg or cluster_cfg[k_req] is None:
            raise ValueError(
                f"FATAL: Missing required clustering parameter 'md_simulation.clustering.{k_req}' in md.yaml!"
            )

    if algo == "hierarchical" and float(cluster_cfg["rmsd_cutoff"]) <= 0.0:
        raise ValueError(
            f"FATAL: 'rmsd_cutoff' must be positive for hierarchical clustering, got: {cluster_cfg['rmsd_cutoff']}"
        )

    if algo in ["kmeans", "auto_kmeans"] and int(cluster_cfg["n_clusters"]) < 1:
        raise ValueError(
            f"FATAL: 'n_clusters' must be a positive integer >= 1, got: {cluster_cfg['n_clusters']}"
        )

    if (
        float(cluster_cfg["variance_cutoff"]) <= 0.0
        or float(cluster_cfg["variance_cutoff"]) > 1.0
    ):
        raise ValueError(
            f"FATAL: 'variance_cutoff' must be a fraction within (0.0, 1.0], got: {cluster_cfg['variance_cutoff']}"
        )

    if float(cluster_cfg["temp_k"]) <= 0.0:
        raise ValueError(
            f"FATAL: 'temp_k' must be a positive absolute temperature in Kelvin, got: {cluster_cfg['temp_k']}"
        )

    plateau_cfg = cluster_cfg.get("plateau_detection")
    if not plateau_cfg or not isinstance(plateau_cfg, dict):
        raise ValueError(
            "FATAL CONFIGURATION ERROR: Missing required configuration section "
            "'md_simulation.clustering.plateau_detection' in configs/md.yaml!\n"
            "Fallbacks and default values are strictly forbidden by Fail-Fast law. "
            "Declare 'plateau_detection' explicitly in configs/md.yaml under 'clustering'."
        )

    required_plateau_keys = [
        "enabled",
        "paradigm",
        "smoothing_tau_ns",
        "significance_z_crit",
        "dilation_buffer_ns",
        "min_dwell_ns",
        "enforce_cluster_continuity",
    ]
    for p_key in required_plateau_keys:
        if p_key not in plateau_cfg or plateau_cfg[p_key] is None or plateau_cfg[p_key] == "":
            raise ValueError(
                f"FATAL CONFIGURATION ERROR: Missing required parameter 'md_simulation.clustering.plateau_detection.{p_key}' "
                "in configs/md.yaml! All parameters must be explicitly specified."
            )

    p_paradigm = str(plateau_cfg["paradigm"]).strip().lower()
    valid_paradigms = ["adaptive_zscore", "macro_step", "running_range", "inverted_vel"]
    if p_paradigm not in valid_paradigms:
        raise ValueError(
            f"FATAL: Invalid plateau_detection.paradigm: '{p_paradigm}'. Must be one of {valid_paradigms}."
        )

    if float(plateau_cfg["smoothing_tau_ns"]) <= 0.0:
        raise ValueError("FATAL: plateau_detection.smoothing_tau_ns must be positive!")
    if float(plateau_cfg["significance_z_crit"]) <= 0.0:
        raise ValueError("FATAL: plateau_detection.significance_z_crit must be positive!")
    if float(plateau_cfg["dilation_buffer_ns"]) < 0.0:
        raise ValueError("FATAL: plateau_detection.dilation_buffer_ns must be non-negative!")
    if float(plateau_cfg["min_dwell_ns"]) <= 0.0:
        raise ValueError("FATAL: plateau_detection.min_dwell_ns must be positive!")

    return cluster_cfg


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
    if not method or method not in ["resp", "bcc", "ashgc"]:
        raise ValueError(
            f"FATAL: Invalid or missing 'md_simulation.ligand_parametrization.method': '{method}'. "
            "Supported methods are 'ashgc' (OpenFF GNN), 'resp' (ORCA 6 quantum ESP), or 'bcc' (AM1-BCC)."
        )

    ff = param_cfg.get("forcefield")
    valid_forcefields = [
        "gaff",
        "gaff2",
        "openff-2.3.0",
        "openff-2.2.1",
        "openff-2.1.0",
    ]
    if not ff or ff not in valid_forcefields:
        raise ValueError(
            f"FATAL: Invalid or missing 'md_simulation.ligand_parametrization.forcefield': '{ff}'. "
            f"Supported small molecule force fields are {valid_forcefields}."
        )

    if method == "ashgc" and not ff.startswith("openff"):
        raise ValueError(
            f"FATAL: Charge method 'ashgc' (OpenFF GNN) requires an OpenFF force field, but got '{ff}'!"
        )

    if method == "resp":
        orca_cfg = param_cfg.get("orca")
        if not orca_cfg or not isinstance(orca_cfg, dict):
            raise ValueError(
                "FATAL: 'md_simulation.ligand_parametrization.orca' section is required when method='resp'!"
            )
        required_orca_keys = ["nprocs", "maxcore_mb"]
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
        nprocs = int(orca["nprocs"])
        maxcore = int(orca["maxcore_mb"])
        total_mem_mb = (nprocs * maxcore) + 2048
        args.extend(
            [
                "--orca_nprocs",
                str(nprocs),
                "--orca_maxcore",
                str(maxcore),
                "--orca_total_mem_mb",
                str(total_mem_mb),
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
