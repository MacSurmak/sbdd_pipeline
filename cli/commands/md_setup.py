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

    # 2. Check cache: identify microstates lacking GAFF2 / AM1-BCC topology
    missing_states = []
    for state_id, state in unique_states.items():
        mol2_ok = state.path_mol2 and Path(state.path_mol2).exists() and Path(state.path_mol2).stat().st_size > 0
        frcmod_ok = state.path_frcmod and Path(state.path_frcmod).exists() and Path(state.path_frcmod).stat().st_size > 0
        if not (mol2_ok and frcmod_ok):
            missing_states.append(state)

    # 3. Lazy On-Demand Parametrization: dispatch quantum calculations strictly for missing microstates
    if missing_states:
        console.print(
            f"[bold cyan]Lazy Parametrization: {len(missing_states)} chemical state(s) require GAFF2/AM1-BCC topology...[/bold cyan]"
        )
        param_plan = []
        from rdkit import Chem
        from rdkit.Chem import AllChem

        for state in missing_states:
            # Prioritize candidate 3D conformers that already contain valid physical hydrogens without zero coordinates
            cand_sdfs = []
            if state.poses:
                for p in state.poses:
                    p_p = p.path_sdf or p.path
                    if p_p and Path(p_p).exists() and Path(p_p).stat().st_size > 0:
                        cand_sdfs.append(Path(p_p).resolve())
            if state.path_sdf and Path(state.path_sdf).exists() and Path(state.path_sdf).stat().st_size > 0:
                cand_sdfs.append(Path(state.path_sdf).resolve())

            best_sdf = None
            for cand in cand_sdfs:
                suppl = Chem.SDMolSupplier(str(cand), removeHs=False, sanitize=False)
                m = next((mol for mol in suppl if mol is not None), None)
                if m and any(a.GetAtomicNum() == 1 for a in m.GetAtoms()):
                    conf = m.GetConformer()
                    has_zero = any(
                        abs(conf.GetAtomPosition(i).x) < 1e-4
                        and abs(conf.GetAtomPosition(i).y) < 1e-4
                        and abs(conf.GetAtomPosition(i).z) < 1e-4
                        for i in range(m.GetNumAtoms())
                    )
                    if not has_zero:
                        best_sdf = cand
                        break

            if best_sdf is None:
                if not cand_sdfs:
                    session.close()
                    raise FileNotFoundError(
                        f"FATAL: SubstanceState '{state.state_hash}' has no reference SDF deliverables."
                    )
                sub_sdf_resolved = cand_sdfs[0]
                suppl = Chem.SDMolSupplier(str(sub_sdf_resolved), removeHs=False, sanitize=False)
                m_base = next((mol for mol in suppl if mol is not None), None)
                if m_base is None:
                    raise ValueError(f"FATAL: Could not read molecule from {sub_sdf_resolved}")

                m_base = Chem.AddHs(m_base, addCoords=True)
                mp = AllChem.MMFFGetMoleculeProperties(m_base)
                if mp:
                    ff = AllChem.MMFFGetMoleculeForceField(m_base, mp)
                    for a in m_base.GetAtoms():
                        if a.GetAtomicNum() != 1:
                            ff.AddFixedPoint(a.GetIdx())
                    ff.Minimize(maxIts=500)
                else:
                    ff = AllChem.UFFGetMoleculeForceField(m_base)
                    if ff:
                        for a in m_base.GetAtoms():
                            if a.GetAtomicNum() != 1:
                                ff.AddFixedPoint(a.GetIdx())
                        ff.Minimize(maxIts=500)

                with Chem.SDWriter(str(sub_sdf_resolved)) as w:
                    w.write(m_base)
                best_sdf = sub_sdf_resolved

            # Atomically synchronize canonical state SDF with validated clean 3D structure
            state.path_sdf = str(best_sdf)

            param_plan.append(
                {
                    "name": f"{state.substance.name}_{state.state_hash}",
                    "res_name": state.substance.res_name,
                    "inchikey": state.substance.inchikey,
                    "state_hash": state.state_hash,
                    "charge": state.formal_charge,
                    "sdf": str(best_sdf),
                }
            )

        param_plan_file = cwd / "param_plan.json"
        param_plan_file.write_text(json.dumps(param_plan, indent=2), encoding="utf-8")
        (cwd / "substance_param_receipt.tsv").unlink(missing_ok=True)

        nf_param_args = [
            "--project_name",
            configs.get("project", {}).get("name", cwd.name),
            "--outdir",
            str(cwd.resolve()),
            "--param_plan",
            str(param_plan_file.resolve()),
        ]
        if urgent:
            nf_param_args.extend(["--queue", "urgent"])

        profile_name = "local" if local else "slurm"
        ret = execute_nextflow(
            "parametrize",
            nf_param_args,
            cwd,
            profile=profile_name,
            resume=resume,
            configs=configs,
        )
        if ret != 0:
            session.close()
            console.print("[bold red]FATAL: On-demand state parametrization failed.[/bold red]")
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
            state_to_update = session.query(SubstanceState).filter_by(state_hash=r_hash).first()
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
