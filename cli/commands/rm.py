"""
File: cli/commands/rm.py
Description: CLI handler to remove items from the clean 3NF database with cascade deletion.
"""

from rich.console import Console
from cli.core.db import (
    init_db,
    Target,
    ReceptorFrame,
    Substance,
    SubstanceState,
    Pose,
    System,
    Trajectory,
    ClusteringRun,
    TrajectoryCluster,
    MmpbsaRun,
    DockingRun,
    DockedPose,
)

console = Console()


def remove_library_substances() -> None:
    """
    Removes all imported screening library substances (batch_id != 'crystal_extract'),
    their physical directories under substances/<inchikey>, any matrix DockingRuns,
    and any unreferenced raw SubstanceStates on crystal ligands.
    Strictly preserves crystal/medoid receptors, native ligands (TGI), ready ORCA states, and MD Systems.
    """
    import shutil
    from pathlib import Path

    session = init_db()
    cwd = Path.cwd()

    matrix_runs = session.query(DockingRun).filter_by(mode="matrix").all()
    for dr in matrix_runs:
        session.delete(dr)
    session.flush()

    lib_subs = (
        session.query(Substance)
        .filter(Substance.batch_id != "crystal_extract")
        .all()
    )
    removed_subs = 0
    for sub in lib_subs:
        sub_dir = cwd / "substances" / sub.inchikey
        if sub_dir.exists():
            shutil.rmtree(sub_dir)
        session.delete(sub)
        removed_subs += 1
    session.flush()

    # Re-link any DockedPose on crystal substances that pointed to a raw state back to its canonical full-H state
    from cli.core.chem import get_microstate_identity
    from rdkit import Chem

    for dp in session.query(DockedPose).all():
        st = dp.substance_state
        if st and st.state == "raw":
            mol = Chem.MolFromMolBlock(dp.mol_block, removeHs=False, sanitize=False)
            if mol is not None:
                _, _, _, canon_hash, _, _ = get_microstate_identity(mol)
                ready_st = (
                    session.query(SubstanceState)
                    .filter_by(substance_id=dp.substance_id, state_hash=canon_hash)
                    .first()
                )
                if ready_st and ready_st.id != st.id:
                    dp.substance_state_id = ready_st.id
    session.flush()

    # Prune any unreferenced raw SubstanceStates on remaining crystal substances
    pruned_states = 0
    for st in session.query(SubstanceState).filter_by(state="raw").all():
        if len(st.poses) == 0 and len(st.docked_poses) == 0:
            session.delete(st)
            pruned_states += 1

    session.commit()
    session.close()
    console.print(
        f"[bold green]✓ Removed {removed_subs} library substance(s), {len(matrix_runs)} matrix DockingRun(s), "
        f"and {pruned_states} orphaned raw microstate(s). Crystal systems and ORCA caches preserved![/bold green]"
    )


def remove_entity(name: str) -> None:
    session = init_db()
    deleted = False

    # Check across models in hierarchical order
    target = session.query(Target).filter_by(name=name).first()
    if target:
        session.delete(target)
        console.print(
            f"[bold green]✓ Deleted Target '{name}' and all associated frames/systems (CASCADE).[/bold green]"
        )
        deleted = True

    if not deleted:
        sys_obj = session.query(System).filter_by(name=name).first()
        if sys_obj:
            session.delete(sys_obj)
            console.print(
                f"[bold green]✓ Deleted System '{name}' and all associated trajectories.[/bold green]"
            )
            deleted = True

    if not deleted:
        frame = session.query(ReceptorFrame).filter_by(name=name).first()
        if frame:
            session.delete(frame)
            console.print(
                f"[bold green]✓ Deleted ReceptorFrame '{name}' and associated poses/systems.[/bold green]"
            )
            deleted = True

    if not deleted:
        pose = session.query(Pose).filter_by(name=name).first()
        if pose:
            session.delete(pose)
            console.print(f"[bold green]✓ Deleted Pose '{name}'.[/bold green]")
            deleted = True

    if not deleted:
        traj = session.query(Trajectory).filter_by(name=name).first()
        if traj:
            session.delete(traj)
            console.print(f"[bold green]✓ Deleted Trajectory '{name}'.[/bold green]")
            deleted = True

    if not deleted:
        c_run = session.query(ClusteringRun).filter_by(name=name).first()
        if c_run:
            session.delete(c_run)
            console.print(
                f"[bold green]✓ Deleted ClusteringRun '{name}' and associated clusters (CASCADE).[/bold green]"
            )
            deleted = True

    if not deleted:
        d_run = session.query(DockingRun).filter_by(name=name).first()
        if d_run:
            session.delete(d_run)
            console.print(
                f"[bold green]✓ Deleted DockingRun '{name}' and all associated docked modes (CASCADE).[/bold green]"
            )
            deleted = True

    if not deleted:
        st = session.query(SubstanceState).filter_by(state_hash=name).first()
        if st:
            session.delete(st)
            console.print(
                f"[bold green]✓ Deleted SubstanceState '{name}'.[/bold green]"
            )
            deleted = True

    if not deleted:
        sub = (
            session.query(Substance)
            .filter((Substance.name == name) | (Substance.inchikey == name))
            .first()
        )
        if sub:
            session.delete(sub)
            console.print(
                f"[bold green]✓ Deleted Chemical Substance '{name}'.[/bold green]"
            )
            deleted = True

    if deleted:
        session.commit()
    else:
        console.print(f"[bold red]✗ Item '{name}' not found in database.[/bold red]")

    session.close()
