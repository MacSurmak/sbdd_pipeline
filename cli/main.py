"""
File: cli/main.py
Description: Main CLI entrypoint for `sbdd` orchestrator tool.
"""

import argparse
import sys
import shutil
import subprocess
from pathlib import Path
from rich.console import Console
from rich.table import Table

console = Console()


def main():
    parser = argparse.ArgumentParser(
        prog="sbdd", description="SBDD Pipeline Nextflow HPC Orchestrator CLI"
    )
    subparsers = parser.add_subparsers(dest="subcommand", help="Available subcommands")

    # 0. add (Database Registration)
    p_add = subparsers.add_parser(
        "add", help="Register external files into the project database"
    )
    add_subs = p_add.add_subparsers(dest="add_type", required=True)

    # 0.1 add structure / receptor
    p_add_rec = add_subs.add_parser(
        "structure",
        aliases=["receptor", "complex"],
        help="Import macromolecular structure(s) into a Target",
    )
    p_add_rec.add_argument(
        "path",
        nargs="+",
        help="Path(s) or glob to structure file(s) (.pdb, .cif, .mae)",
    )
    p_add_rec.add_argument(
        "--target", help="Parent target name (e.g. RyR1, TSPO). Defaults to file stem."
    )
    p_add_rec.add_argument("--name", help="Unique identifier for this structural frame")
    p_add_rec.add_argument(
        "--cofactors",
        nargs="*",
        default=[],
        help="Cofactor residue names to retain (e.g. ATP, NADH, FAD)",
    )
    p_add_rec.add_argument(
        "--ions",
        nargs="*",
        default=[],
        help="Explicit ion residue names to retain (e.g. MG, ZN, CA, CL)",
    )
    p_add_rec.add_argument(
        "--ligand", help="Native ligand residue name inside structure (e.g. KVR, BTN)"
    )
    p_add_rec.add_argument(
        "--fasta", help="Explicit path to reference FASTA sequence for the target"
    )

    # 0.2 add fasta
    p_add_fasta = add_subs.add_parser(
        "fasta", help="Import reference FASTA sequence for a Target"
    )
    p_add_fasta.add_argument("path", help="Path to the .fasta file")

    # 0.3 add ligand
    p_add_lig = add_subs.add_parser(
        "ligand", help="Register a ligand library or individual molecule"
    )
    p_add_lig.add_argument("path", help="Path to the ligand file (.sdf)")
    p_add_lig.add_argument(
        "--name", help="Identifier for ligand (optional for multi-molecule SDF)"
    )
    p_add_lig.add_argument(
        "--force", action="store_true", help="Overwrite if name exists"
    )

    # 0.4 add traj
    p_add_traj = add_subs.add_parser(
        "traj", help="Register an existing MD trajectory replicate linked to a System"
    )
    p_add_traj.add_argument(
        "--system",
        required=True,
        help="Name of the System in DB (e.g. 9ZKE_crystal_holo)",
    )
    p_add_traj.add_argument(
        "--replicate",
        type=int,
        default=1,
        help="Replicate index (1, 2, 3...)",
    )
    p_add_traj.add_argument(
        "--raw-dir",
        help="Path to directory containing raw chunk files md_*.xtc and sys_ref.tpr",
    )
    p_add_traj.add_argument(
        "--xtc",
        help="Path to raw .xtc file or primary chunk (optional if --raw-dir is provided)",
    )
    p_add_traj.add_argument(
        "--tpr",
        required=True,
        help="Path to .tpr file (e.g. sys_ref.tpr)",
    )
    p_add_traj.add_argument(
        "--length_ns",
        type=float,
        default=0.0,
        help="Total simulation length in nanoseconds (e.g. 100.0)",
    )
    p_add_traj.add_argument(
        "--name",
        help="Unique name for this trajectory (default: <system>_rep<replicate>)",
    )
    p_add_traj.add_argument(
        "--clean-xtc",
        help="Explicit path to pre-cleaned .xtc (optional, if already cleaned)",
    )
    p_add_traj.add_argument(
        "--ref-pdb",
        help="Explicit path to dry reference .pdb (optional, if already cleaned)",
    )
    p_add_traj.add_argument("--force", action="store_true", help="Overwrite if exists")
    # 0.4 tree
    p_tree = subparsers.add_parser(
        "tree", aliases=["ls"], help="Show database lineage tree"
    )

    # 0.4 rm
    p_rm = subparsers.add_parser("rm", help="Remove an item from the database by name")
    p_rm.add_argument("target", help="Name of the item to remove")

    # Common options for orchestration
    def add_common_wf_args(p):
        p.add_argument("--local", action="store_true", help="Run locally without Slurm")
        p.add_argument(
            "--no-resume",
            action="store_true",
            help="Start from scratch without -resume",
        )
        p.add_argument(
            "--clean-scratch",
            action="store_true",
            help="Delete /shared/work upon completion (default: False)",
        )
        p.add_argument(
            "--urgent",
            action="store_true",
            help="Route tasks to urgent Slurm partition",
        )

    # 1. init
    p_init = subparsers.add_parser("init", help="Initialize new project workspace")
    p_init.add_argument(
        "directory",
        nargs="?",
        default=".",
        help="Target project directory (default: current)",
    )
    p_init.add_argument(
        "--name", help="Explicit project name (default: auto-generated random name)"
    )
    p_init.add_argument(
        "--force", action="store_true", help="Overwrite existing templates"
    )

    # 1.5 heal
    p_heal = subparsers.add_parser(
        "heal", help="Physically repair missing loops in broken crystal frames via AI"
    )
    p_heal.add_argument("--target", help="Limit healing to specific Target")
    add_common_wf_args(p_heal)

    # 2. prep / prepare
    p_prep = subparsers.add_parser(
        "prep",
        aliases=["prepare"],
        help="Prepare raw receptor frames and in-situ poses into ready systems",
    )
    p_prep.add_argument(
        "--target",
        help="Prepare all raw frames for a specific Target",
    )
    p_prep.add_argument(
        "--frames",
        help="Comma-separated list of specific ReceptorFrame names to prepare",
    )
    add_common_wf_args(p_prep)

    # 3. dock / screen
    p_dock = subparsers.add_parser(
        "dock", aliases=["screen"], help="Run distributed docking into ReceptorFrames"
    )
    p_dock.add_argument(
        "--target", help="Parent target name to dock into its ready frames"
    )
    p_dock.add_argument(
        "--frames",
        "--receptors",
        dest="frames",
        help="Comma-separated list of specific ready ReceptorFrame names",
    )
    p_dock.add_argument(
        "--ligands",
        help="Comma-separated list of ligand names from DB (optional in --self / --redock mode)",
    )
    p_dock.add_argument(
        "--self",
        "--redock",
        dest="self_dock",
        action="store_true",
        help="Run 1-to-1 self-docking (redocking) of each receptor with its own native ligand without cross-docking matrix",
    )
    p_dock.add_argument("--name", help="Name for the output PoseSet in DB")
    add_common_wf_args(p_dock)

    # 3.1 select-hits
    p_hits = subparsers.add_parser(
        "select-hits",
        help="Select top candidate hits from docking rankings and promote to Systems",
    )
    p_hits.add_argument(
        "--poses", help="Name of PoseSet in DB to select hits from (default: latest)"
    )
    p_hits.add_argument(
        "--top", type=int, default=5, help="Number of top hits to extract (default: 5)"
    )
    p_hits.add_argument(
        "--system-prefix",
        default="hit",
        help="Prefix for generated Holo systems (default: hit)",
    )

    # 4. rescore
    p_rescore = subparsers.add_parser(
        "rescore",
        help="Run standalone Prime MM-GBSA free energy rescoring on docked poses",
    )
    p_rescore.add_argument(
        "--use-whitelist",
        action="store_true",
        help="Filter poses using funnel whitelist prior to MM-GBSA",
    )
    add_common_wf_args(p_rescore)

    # 5. md-init
    p_mdinit = subparsers.add_parser(
        "md-init", help="Initialize editable MD protocol playlist"
    )

    # 5.5 md-setup
    p_mdsetup = subparsers.add_parser(
        "md-setup", help="Prepare MD topologies for raw systems (tLEaP/Declash)"
    )
    p_mdsetup.add_argument("--system", help="Specific System to setup")
    p_mdsetup.add_argument(
        "--target", help="Setup all raw systems for a specific Target"
    )
    p_mdsetup.add_argument(
        "--force",
        action="store_true",
        help="Force re-parametrization of chemical states using the configured method (e.g. ORCA RESP)",
    )
    add_common_wf_args(p_mdsetup)

    # 6. md / dynamics
    p_dyn = subparsers.add_parser(
        "dynamics",
        aliases=["md", "md-run"],
        help="Execute chunked MD simulations for a ready System using protocol.txt playlist",
    )
    p_dyn.add_argument(
        "--system", help="Name of the physical System to simulate (from 'sbdd tree')"
    )
    p_dyn.add_argument(
        "--target", help="Simulate all ready Systems belonging to a specific Target"
    )
    add_common_wf_args(p_dyn)

    # 7. traj-mmpbsa
    p_mmpbsa = subparsers.add_parser(
        "traj-mmpbsa",
        help="Run standalone MM-PBSA free energy calculation on clean trajectories",
    )
    p_mmpbsa.add_argument("--systems", help="Limit to specific System names")
    p_mmpbsa.add_argument("--target", help="Limit to specific Target")
    p_mmpbsa.add_argument(
        "--traj", help="Explicit comma-separated list of Trajectory names"
    )
    add_common_wf_args(p_mmpbsa)

    # 7. md-analyze
    p_mda = subparsers.add_parser(
        "md-analyze", help="Run trajectory PBC fix, metrics, and MM-PBSA"
    )
    p_mda.add_argument("--traj", help="Direct path to custom .xtc file")
    p_mda.add_argument("--tpr", help="Direct path to custom .tpr file")
    p_mda.add_argument("--topol", help="Direct path to custom topol.top or PDB")
    p_mda.add_argument("--ndx", help="Direct path to custom index.ndx file")
    add_common_wf_args(p_mda)

    # 7.0 traj-flag
    p_flag = subparsers.add_parser(
        "traj-flag",
        help="Flag or disable a trajectory (e.g. ligand drifted/unbound) to exclude from clustering",
    )
    p_flag.add_argument("trajectory", help="Name of the Trajectory record in DB")
    p_flag.add_argument(
        "--status",
        choices=["cleaned", "drifted", "failed", "disabled", "raw"],
        default="drifted",
        help="Target QC status (default: drifted)",
    )
    p_flag.add_argument(
        "--disable", action="store_true", help="Shortcut to set status='disabled'"
    )
    p_flag.add_argument(
        "--enable", action="store_true", help="Shortcut to restore status='cleaned'"
    )

    # 7.1 traj-clean
    p_clean_tr = subparsers.add_parser(
        "traj-clean",
        help="Remove solvent, fix PBC, center on protein, and fit trajectory for all uncleaned runs",
    )
    p_clean_tr.add_argument(
        "--systems",
        help="Limit to specific System names",
    )
    p_clean_tr.add_argument(
        "--target",
        help="Limit to specific Target",
    )
    p_clean_tr.add_argument(
        "--traj", help="Explicit comma-separated list of Trajectory names"
    )
    add_common_wf_args(p_clean_tr)

    # 7.2 traj-metrics
    p_metrics = subparsers.add_parser(
        "traj-metrics",
        help="Calculate RMSD (protein, ligand, and cofactors), RMSF, and contacts for cleaned trajectories",
    )
    p_metrics.add_argument("--systems", help="Limit to specific System names")
    p_metrics.add_argument("--target", help="Limit to specific Target")
    p_metrics.add_argument(
        "--traj", help="Explicit comma-separated list of Trajectory names"
    )
    add_common_wf_args(p_metrics)

    # 8. traj-cluster
    p_cluster = subparsers.add_parser(
        "traj-cluster",
        help="Cluster clean MD trajectory ensembles using PCA/TICA + UMAP macro-basins and register centroids in 3NF DB",
    )
    p_cluster.add_argument(
        "--systems",
        help="Limit clustering to specific System names",
    )
    p_cluster.add_argument(
        "--target",
        help="Limit clustering to a specific Target",
    )
    p_cluster.add_argument(
        "--traj", help="Explicit comma-separated list of Trajectory names"
    )
    p_cluster.add_argument(
        "--align-group",
        help="Basis group for coordinate superimposition (default from md.yaml, e.g. Backbone)",
    )
    p_cluster.add_argument(
        "--cluster-group",
        "--group",
        dest="cluster_group",
        help="Feature group for manifold reduction (default from md.yaml, e.g. Pocket, Backbone, Ligand)",
    )
    p_cluster.add_argument(
        "--dim-reduction",
        choices=["pca", "tica"],
        help="Dimensionality reduction engine: pca or tica (default from md.yaml)",
    )
    p_cluster.add_argument(
        "--n-clusters",
        "-k",
        type=int,
        help="Number of physical macro-basins to partition (default from md.yaml)",
    )
    p_cluster.add_argument(
        "--min-population",
        type=float,
        help="Population threshold (%) for automated CLI recommendation (default from md.yaml)",
    )
    p_cluster.add_argument(
        "--auto",
        "--yes",
        dest="auto",
        action="store_true",
        help="Headless non-interactive execution: automatically promote recommended clusters into DB",
    )
    p_cluster.add_argument(
        "--force",
        action="store_true",
        help="Force re-clustering and bypass Nextflow caching",
    )
    add_common_wf_args(p_cluster)

    # 8.1 traj-cluster-promote (Dedicated, non-blocking cluster promotion command)
    p_promote = subparsers.add_parser(
        "traj-cluster-promote",
        help="Promote specific cluster centroids from an existing manifest into SQLite 3NF without re-running MDAnalysis",
    )
    p_promote.add_argument("--manifest", required=True, help="Path to clusters_manifest.json")
    p_promote.add_argument("--ranks", required=True, help="Comma-separated ranks to promote (e.g. 1, 3)")

    # 9. clean
    p_clean = subparsers.add_parser(
        "clean", help="Clean /shared/work scratch directory for current project"
    )

    # 9. status
    p_stat = subparsers.add_parser(
        "status", help="Inspect Slurm queue and Nextflow job status"
    )

    # 10. watch / progress
    p_watch = subparsers.add_parser(
        "watch",
        aliases=["progress", "top"],
        help="Monitor real-time task progress via cluster UDP telemetry",
    )
    p_watch.add_argument(
        "--port",
        type=int,
        default=8088,
        help="Telemetry relay UDP port (default: 8088)",
    )

    args = parser.parse_args()

    if not args.subcommand:
        parser.print_help()
        sys.exit(0)

    subcmd = args.subcommand

    if subcmd == "add":
        from cli.commands.add import add_receptor, add_ligand, add_traj, add_fasta

        if args.add_type in ["structure", "receptor", "complex"]:
            add_receptor(
                paths=args.path if isinstance(args.path, list) else [args.path],
                target=args.target,
                name=args.name,
                cofactors=args.cofactors,
                ions=args.ions,
                ligand=args.ligand,
                fasta=args.fasta,
            )
        elif args.add_type == "fasta":
            add_fasta(path=args.path)
        elif args.add_type == "ligand":
            add_ligand(path=args.path, name=args.name, force=args.force)
        elif args.add_type == "traj":
            add_traj(
                system=args.system,
                tpr=args.tpr,
                xtc=args.xtc,
                raw_dir=args.raw_dir,
                replicate=args.replicate,
                name=args.name,
                length_ns=args.length_ns,
                clean_xtc=args.clean_xtc,
                ref_pdb=args.ref_pdb,
            )
    elif subcmd in ["tree", "ls"]:
        from cli.commands.tree import print_tree

        print_tree()
    elif subcmd == "rm":
        from cli.commands.rm import remove_entity

        remove_entity(name=args.target)
    elif subcmd == "init":
        from cli.commands.project import run_init

        run_init(directory=args.directory, name=args.name, force=args.force)
    elif subcmd == "heal":
        from cli.commands.heal import run_heal

        run_heal(
            target=args.target,
            urgent=args.urgent,
            local=args.local,
            resume=not args.no_resume,
        )
    elif subcmd == "md-init":
        from cli.commands.project import run_md_init

        run_md_init()
    elif subcmd == "status":
        from cli.commands.project import run_status

        run_status()
    elif subcmd == "clean":
        from cli.commands.project import run_clean

        run_clean()
    elif subcmd in ["prep", "prepare"]:
        from cli.commands.prep import run_prep

        run_prep(
            target=args.target,
            frames=(
                [f.strip() for f in args.frames.split(",") if f.strip()]
                if args.frames
                else None
            ),
            urgent=args.urgent,
            local=args.local,
            resume=not args.no_resume,
        )
    elif subcmd in ["dock", "screen"]:
        from cli.commands.dock import run_dock

        run_dock(
            target=args.target,
            frame_names=(
                [f.strip() for f in args.frames.split(",") if f.strip()]
                if args.frames
                else None
            ),
            ligand_names=(
                [n.strip() for n in args.ligands.split(",")] if args.ligands else None
            ),
            self_dock=args.self_dock,
            name=args.name,
            urgent=args.urgent,
            local=args.local,
            resume=not args.no_resume,
        )
    elif subcmd == "rescore":
        from cli.commands.rescore import run_rescore

        run_rescore(
            use_whitelist=args.use_whitelist,
            urgent=args.urgent,
            local=args.local,
            resume=not args.no_resume,
        )
    elif subcmd == "select-hits":
        from cli.commands.select_hits import run_select_hits

        run_select_hits(
            poses_name=args.poses,
            top_n=args.top,
            system_prefix=args.system_prefix,
        )
    elif subcmd == "md-setup":
        from cli.commands.md_setup import run_md_setup

        run_md_setup(
            system=args.system,
            target=args.target,
            force=args.force,
            urgent=args.urgent,
            local=args.local,
            resume=not args.no_resume,
        )
    elif subcmd in ["dynamics", "md", "md-run"]:
        from cli.commands.dynamics import run_dynamics

        run_dynamics(
            system=args.system,
            target=args.target,
            urgent=args.urgent,
            local=args.local,
            resume=not args.no_resume,
        )
    elif subcmd == "md-analyze":
        from cli.commands.traj import run_md_analyze

        run_md_analyze(
            traj=args.traj,
            tpr=args.tpr,
            topol=args.topol,
            ndx=args.ndx,
            urgent=args.urgent,
            local=args.local,
            resume=not args.no_resume,
            clean_scratch=args.clean_scratch,
        )
    elif subcmd == "traj-flag":
        from cli.commands.traj import run_flag_traj

        run_flag_traj(
            trajectory=args.trajectory,
            status=args.status,
            disable=args.disable,
            enable=args.enable,
        )
    elif subcmd == "traj-clean":
        from cli.commands.traj import run_clean_traj

        run_clean_traj(
            systems=(
                [s.strip() for s in args.systems.split(",") if s.strip()]
                if args.systems
                else None
            ),
            target=args.target,
            traj_names=(
                [t.strip() for t in args.traj.split(",") if t.strip()]
                if args.traj
                else None
            ),
            urgent=args.urgent,
            local=args.local,
            resume=not args.no_resume,
        )
    elif subcmd == "traj-mmpbsa":
        from cli.commands.traj import run_mmpbsa

        run_mmpbsa(
            systems=(
                [s.strip() for s in args.systems.split(",") if s.strip()]
                if args.systems
                else None
            ),
            target=args.target,
            traj_names=(
                [t.strip() for t in args.traj.split(",") if t.strip()]
                if args.traj
                else None
            ),
            urgent=args.urgent,
            local=args.local,
            resume=not args.no_resume,
        )
    elif subcmd == "traj-metrics":
        from cli.commands.traj import run_metrics

        run_metrics(
            systems=(
                [s.strip() for s in args.systems.split(",") if s.strip()]
                if args.systems
                else None
            ),
            target=args.target,
            traj_names=(
                [t.strip() for t in args.traj.split(",") if t.strip()]
                if args.traj
                else None
            ),
            urgent=args.urgent,
            local=args.local,
            resume=not args.no_resume,
        )
    elif subcmd == "traj-cluster":
        from cli.commands.traj import run_cluster

        run_cluster(
            systems=(
                [s.strip() for s in args.systems.split(",") if s.strip()]
                if args.systems
                else None
            ),
            target=args.target,
            traj_names=(
                [t.strip() for t in args.traj.split(",") if t.strip()]
                if args.traj
                else None
            ),
            align_group=args.align_group,
            cluster_group=args.cluster_group,
            dim_reduction=args.dim_reduction,
            n_clusters=args.n_clusters,
            min_population=args.min_population,
            auto=args.auto,
            force=args.force,
            urgent=args.urgent,
            local=args.local,
            resume=not args.no_resume,
        )
    elif subcmd == "traj-cluster-promote":
        from cli.commands.traj import run_cluster_promote

        run_cluster_promote(
            manifest_path=args.manifest,
            ranks=[int(x.strip().replace("#", "")) for x in args.ranks.split(",") if x.strip()],
        )
    elif subcmd in ["watch", "progress", "top"]:
        from cli.core.monitor import run_monitor

        run_monitor(port=args.port)


if __name__ == "__main__":
    main()
