"""
File: cli/commands/md_setup.py
Description: CLI handler for 'sbdd md-setup'.
Gathers 'raw' Systems, generates topologies, and atomically registers them as 'ready'.
"""

import json
import sys
from pathlib import Path
from rich.console import Console
from cli.core.db import init_db, Target, System, Substance, SubstanceState, Pose
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def run_md_setup(
    system: str | None = None,
    target: str | None = None,
    force: bool = False,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
) -> None:
    cwd = Path.cwd()
    session = init_db()
    configs = load_project_configs(cwd)

    query = session.query(System).filter_by(state="raw")
    if system:
        query = query.filter_by(name=system)
    elif target:
        target_obj = session.query(Target).filter_by(name=target).first()
        if not target_obj:
            console.print(f"[bold red]FATAL: Target '{target}' not found.[/bold red]")
            sys.exit(1)
        query = query.filter_by(target_id=target_obj.id)

    raw_systems = query.all()
    if not raw_systems:
        console.print(
            "[green]No raw Systems found to setup. Topologies are up to date![/green]"
        )
        session.close()
        return

    # 1. Collect all unique non-ion chemical microstates required across candidate systems
    unique_states = {}
    for s in raw_systems:
        if s.pose and s.pose.substance.substance_type != "ion":
            state = s.pose.substance_state
            if state:
                unique_states[state.id] = state
        for c in s.cofactors:
            if c.substance.substance_type != "ion":
                state = c.substance_state
                if state:
                    unique_states[state.id] = state

    # 2. Check cache: identify microstates lacking valid topology matching target forcefield & method
    from cli.core.config_loader import validate_ligand_parametrization_config

    param_cfg = validate_ligand_parametrization_config(configs)
    target_ff = param_cfg["forcefield"].lower()
    target_method = param_cfg["method"].lower()

    missing_states = []
    for state_id, state in unique_states.items():
        if force:
            missing_states.append(state)
            continue
        mol2_ok = (
            state.path_mol2
            and Path(state.path_mol2).exists()
            and Path(state.path_mol2).stat().st_size > 0
        )
        frcmod_ok = (
            state.path_frcmod
            and Path(state.path_frcmod).exists()
            and Path(state.path_frcmod).stat().st_size > 0
        )
        cache_valid = False
        if mol2_ok and frcmod_ok:
            frcmod_content = Path(state.path_frcmod).read_text(encoding="utf-8", errors="ignore")
            is_openff_cached = "OpenFF" in frcmod_content
            if target_ff.startswith("openff"):
                ff_match = target_ff in frcmod_content.lower()
                method_match = f"method: {target_method}" in frcmod_content.lower()
                cache_valid = is_openff_cached and ff_match and method_match
            else:
                cache_valid = not is_openff_cached

        if not cache_valid:
            missing_states.append(state)

    # 3. Lazy On-Demand Parametrization: dispatch quantum calculations strictly for missing microstates
    if missing_states:
        method_label = "ORCA RESP" if param_cfg["method"] == "resp" else ("AshGC GNN" if param_cfg["method"] == "ashgc" else "AM1-BCC")
        console.print(
            f"[bold cyan]Lazy Parametrization ({method_label}): {len(missing_states)} chemical state(s) require {param_cfg['forcefield'].upper()} topology...[/bold cyan]"
        )
        param_plan = []
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from cli.core.db import system_cofactors

        # Globally resolve deterministic cofactor / secondary component slots from 3NF database
        cofactor_states_query = (
            session.query(SubstanceState.id)
            .join(Pose, Pose.substance_state_id == SubstanceState.id)
            .join(system_cofactors, system_cofactors.c.pose_id == Pose.id)
            .order_by(SubstanceState.id)
            .distinct()
            .all()
        )
        cofactor_order = [r[0] for r in cofactor_states_query]
        # Ahead-of-time immune prefix slots: active ligand gets '1' (cap: 310), cofactors get 2-prefix slots (cap: 124)
        cofactor_prefixes = ["6", "8", "0", "Q", "U", "Y"]

        for state in missing_states:
            is_cofactor = state.id in cofactor_order
            if is_cofactor:
                c_idx = cofactor_order.index(state.id)
                if c_idx >= len(cofactor_prefixes):
                    raise OverflowError(
                        f"FATAL: Exceeded maximum concurrent cofactors ({len(cofactor_prefixes)}) supported by disjoint AMBER namespaces!"
                    )
                assigned_prefix = cofactor_prefixes[c_idx]
            else:
                assigned_prefix = "1"
            # Strict Fail-Fast: reference 3D conformer must exist with physical hydrogens from PrepWizard/LigPrep
            if (
                not state.path_sdf
                or not Path(state.path_sdf).exists()
                or Path(state.path_sdf).stat().st_size == 0
            ):
                session.close()
                raise FileNotFoundError(
                    f"FATAL: SubstanceState '{state.state_hash}' has no valid prepared reference SDF: {state.path_sdf}. "
                    "Molecule must be prepared via PrepWizard or LigPrep before topology setup."
                )

            state_sdf = Path(state.path_sdf).resolve()
            suppl = Chem.SDMolSupplier(str(state_sdf), removeHs=False, sanitize=False)
            m = next((mol for mol in suppl if mol is not None), None)
            if m is None:
                session.close()
                raise ValueError(
                    f"FATAL: Could not parse chemical structure from {state_sdf}"
                )

            if not any(a.GetAtomicNum() == 1 for a in m.GetAtoms()):
                session.close()
                raise ValueError(
                    f"FATAL: SubstanceState '{state.state_hash}' in {state_sdf} lacks explicit hydrogens! "
                    "All hydrogens must be assigned deterministically by PrepWizard or LigPrep at target pH."
                )

            param_plan.append(
                {
                    "name": f"{state.substance.name}_{state.state_hash}",
                    "res_name": state.substance.res_name,
                    "inchikey": state.substance.inchikey,
                    "state_hash": state.state_hash,
                    "charge": state.formal_charge,
                    "method": param_cfg["method"],
                    "forcefield": param_cfg["forcefield"],
                    "atom_type_prefix": assigned_prefix,
                    "sdf": str(state_sdf),
                }
            )

        param_plan_file = cwd / "param_plan.json"
        param_plan_file.write_text(json.dumps(param_plan, indent=2), encoding="utf-8")
        (cwd / "substance_param_receipt.tsv").unlink(missing_ok=True)

        nf_param_args = build_nextflow_args(configs, cwd, entrypoint="parametrize")
        nf_param_args.extend(["--param_plan", str(param_plan_file.resolve())])
        if urgent:
            nf_param_args.extend(["--queue", "urgent"])

        profile_name = "local" if local else "slurm"
        resume_param = resume and not force
        ret = execute_nextflow(
            "parametrize",
            nf_param_args,
            cwd,
            profile=profile_name,
            resume=resume_param,
            configs=configs,
        )
        if ret != 0:
            session.close()
            console.print(
                "[bold red]FATAL: On-demand state parametrization failed.[/bold red]"
            )
            sys.exit(ret)

        param_receipt = cwd / "substance_param_receipt.tsv"
        if not param_receipt.exists() or param_receipt.stat().st_size == 0:
            session.close()
            raise FileNotFoundError(
                f"FATAL: Expected substance parametrization receipt {param_receipt} not found."
            )

        for line in param_receipt.read_text().splitlines():
            if not line.strip():
                continue
            r_hash, r_mol2, r_frcmod = line.strip().split("\t")
            state_to_update = (
                session.query(SubstanceState).filter_by(state_hash=r_hash).first()
            )
            if state_to_update:
                state_to_update.path_mol2 = str(Path(r_mol2).resolve())
                state_to_update.path_frcmod = str(Path(r_frcmod).resolve())
                state_to_update.state = "ready"

        session.commit()
        console.print(
            f"[bold green]✓ Parametrized {len(missing_states)} chemical state(s). Cached for lifetime of project![/bold green]"
        )

    plan_list = []
    for s in raw_systems:
        frame = s.receptor_frame
        if not frame.path_pdb:
            console.print(
                f"[bold red]FATAL: ReceptorFrame '{frame.name}' missing path_pdb. Run 'sbdd prep' first.[/bold red]"
            )
            sys.exit(1)

        rec_pdb = Path(frame.path_pdb).resolve()
        if not rec_pdb.exists():
            raise FileNotFoundError(f"Prepared receptor PDB missing: {rec_pdb}")

        lig_sdf = None
        lig_resname = None
        ion_sdfs = []
        ion_resnames = []
        cofactor_sdfs = []
        cofactor_resnames = []
        mol2_paths = []
        frcmod_paths = []

        # 1. Active ligand (topology resolved strictly from SubstanceState SSOT)
        if s.pose:
            if s.pose.substance.substance_type != "ion":
                state = s.pose.substance_state
                state_mol2 = state.path_mol2
                state_frcmod = state.path_frcmod
                if not state_mol2 or not state_frcmod or not Path(state_mol2).exists():
                    console.print(
                        f"[bold red]FATAL: SubstanceState '{state.state_hash}' missing GAFF2 topology.[/bold red]"
                    )
                    sys.exit(1)
                mol2_paths.append(str(Path(state_mol2).resolve()))
                frcmod_paths.append(str(Path(state_frcmod).resolve()))

            p_sdf = s.pose.path_sdf or s.pose.path
            if not p_sdf:
                console.print(
                    f"[bold red]FATAL: Pose '{s.pose.name}' missing path_sdf.[/bold red]"
                )
                sys.exit(1)
            lig_sdf = str(Path(p_sdf).resolve())
            lig_resname = s.pose.res_name

        # 2. Strictly separate Ions and Organic Cofactors preserving 3NF entity types
        for c in s.cofactors:
            c_sdf = c.path_sdf or c.path
            if not c_sdf or not Path(c_sdf).exists():
                console.print(
                    f"[bold red]FATAL: System '{s.name}' component '{c.name}' missing SDF deliverable.[/bold red]"
                )
                sys.exit(1)

            sdf_path_resolved = str(Path(c_sdf).resolve())

            if c.substance.substance_type == "ion":
                ion_sdfs.append(sdf_path_resolved)
                ion_resnames.append(c.res_name)
            else:
                cofactor_sdfs.append(sdf_path_resolved)
                cofactor_resnames.append(c.res_name)
                state = c.substance_state
                c_mol2 = state.path_mol2
                c_frcmod = state.path_frcmod
                if not c_mol2 or not c_frcmod or not Path(c_mol2).exists():
                    console.print(
                        f"[bold red]FATAL: Cofactor SubstanceState '{state.state_hash}' missing GAFF2 topology.[/bold red]"
                    )
                    sys.exit(1)
                mol2_paths.append(str(Path(c_mol2).resolve()))
                frcmod_paths.append(str(Path(c_frcmod).resolve()))

        plan_list.append(
            {
                "id": s.name,
                "rec_pdb": str(rec_pdb),
                "lig_sdf": lig_sdf,
                "lig_resname": lig_resname,
                "ion_sdfs": ion_sdfs,
                "ion_resnames": ion_resnames,
                "cofactor_sdfs": cofactor_sdfs,
                "cofactor_resnames": cofactor_resnames,
                "mol2_paths": mol2_paths,
                "frcmod_paths": frcmod_paths,
            }
        )

    plan_file = cwd / "md_setup_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_list, f, indent=2)
    (cwd / "md_setup_receipt.tsv").unlink(missing_ok=True)

    console.print(
        f"[bold cyan]Dispatching topology setup for {len(plan_list)} System(s)...[/bold cyan]"
    )

    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-setup")
    nf_args.extend(["--md_plan", str(plan_file.resolve())])

    if urgent:
        nf_args.extend(["--queue", "urgent"])

    profile_name = "local" if local else "slurm"
    ret_code = execute_nextflow(
        "md-setup", nf_args, cwd, profile=profile_name, resume=resume, configs=configs
    )

    if ret_code != 0:
        console.print("[bold red]Nextflow MD setup failed.[/bold red]")
        sys.exit(ret_code)

    receipt_file = cwd / "md_setup_receipt.tsv"
    if not receipt_file.exists():
        console.print("[bold red]FATAL: md_setup_receipt.tsv not found.[/bold red]")
        sys.exit(1)

    sys_map = {s.name: s for s in raw_systems}
    registered = 0
    lines = receipt_file.read_text().strip().splitlines()
    for line in lines:
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) == 5:
            sys_id, gro, top, ndx, posre_str = parts
            posre_files = posre_str.split(",") if posre_str else []
        elif len(parts) == 4:
            sys_id, gro, top, ndx = parts
            posre_files = []
        else:
            continue

        if sys_id in sys_map:
            s = sys_map[sys_id]
            s.path_gro = gro
            s.path_top = top
            s.path_ndx = ndx

            env = s.environment_dict
            env["posre_files"] = posre_files
            s.environment_dict = env

            s.state = "ready"
            registered += 1

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ MD Setup completed! {registered} System(s) transitioned to 'ready'.[/bold green]"
    )
    console.print(
        "[dim]Topologies compiled. Run 'sbdd md' to launch production dynamics.[/dim]"
    )
