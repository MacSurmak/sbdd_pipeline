"""
File: cli/core/runner.py
Description: Nextflow orchestrator execution wrapper with Slurm monitoring and /shared scratch lifecycle.
"""

import subprocess
import shutil
import sys
from pathlib import Path
from rich.console import Console

console = Console()
FRAMEWORK_DIR = Path(__file__).resolve().parent.parent.parent


def clean_shared_scratch(workdir: str, verbose: bool = True):
    """Safely cleans the project's scratch workspace."""
    scratch_dir = Path(workdir)
    if scratch_dir.exists():
        if verbose:
            console.print(
                f"[yellow]Cleaning temporary files in {scratch_dir}...[/yellow]"
            )
        try:
            shutil.rmtree(scratch_dir)
            if verbose:
                console.print(f"[green]Successfully cleaned {scratch_dir}[/green]")
        except Exception as e:
            console.print(f"[red]Error cleaning {scratch_dir}: {e}[/red]")
    else:
        if verbose:
            console.print(f"[dim]No scratch directory found at {scratch_dir}[/dim]")


def build_unrolled_md_dag(total_chunks: int) -> str:
    includes = [
        "include { GMX_EQUILIBRATION } from '../modules/gromacs/equilibration'",
        "include { TRJCAT            } from '../modules/gromacs/trjcat'",
    ]
    for i in range(1, total_chunks + 1):
        includes.append(
            f"include {{ RUN_MD_CHUNK as RUN_MD_CHUNK_{i} }} from '../modules/gromacs/run_chunk'"
        )

    includes_str = "\n".join(includes)

    code = f"""/* Dynamically generated unrolled DAG for chunked MD */
{includes_str}

workflow RUN_CHUNKED_MD {{
    take:
    ch_system
    ch_md_protocol
    total_ns
    chunk_ns
    replicates

    main:
    GMX_EQUILIBRATION(ch_system, ch_md_protocol)
    
    def num_reps = replicates as int
    ch_state_0 = GMX_EQUILIBRATION.out.equilibrated_state
        .join(ch_system.map{{ m, gro, top, ndx, pos -> [m, top, ndx, pos] }}, by: 0)
        .flatMap {{ meta, eq_gro, eq_cpt, eq_tpr, topol, ndx, posre_files ->
            (1..num_reps).collect {{ rep_id ->
                def rep_meta = meta + [rep: rep_id]
                [ rep_meta, eq_gro, eq_cpt, eq_tpr, topol, ndx, posre_files ]
            }}
        }}

    ch_all_xtcs = Channel.empty()
"""
    for i in range(1, total_chunks + 1):
        code += f"""
    step_{i} = RUN_MD_CHUNK_{i}(ch_state_{i-1}, ch_md_protocol, chunk_ns, {i})
    ch_state_{i} = step_{i}.state.map {{ meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }}
    ch_all_xtcs = ch_all_xtcs.mix(step_{i}.xtc)
"""
    code += f"""
    sys_files = step_1.sys_files

    ch_grouped_xtcs = ch_all_xtcs.groupTuple(by: 0, size: {total_chunks})
    TRJCAT(ch_grouped_xtcs)

    emit:
    raw_trajectory = TRJCAT.out.concat_xtc
    sys_files = sys_files
}}
"""
    return code


def execute_nextflow(
    entrypoint: str,
    extra_args: list,
    project_dir: Path,
    profile: str = "slurm",
    resume: bool = True,
    auto_clean: bool = False,
    configs: dict = None,
) -> int:
    """Invokes Nextflow orchestrator with complete isolation from project directory."""
    import os

    main_nf = FRAMEWORK_DIR / "main.nf"

    project_name = project_dir.name
    if "--project_name" in extra_args:
        idx = extra_args.index("--project_name")
        if idx + 1 < len(extra_args):
            project_name = extra_args[idx + 1]

    # Resolve scratch base safely via config or fallback
    default_scratch = f"/shared/work/{project_name}"
    scratch_base_str = default_scratch
    if configs:
        scratch_base_str = configs.get("project", {}).get("workdir", default_scratch)

    scratch_base = Path(scratch_base_str)
    scratch_nxf_home = scratch_base / ".nextflow"
    scratch_work = scratch_base / "work"
    scratch_nxf_home.mkdir(parents=True, exist_ok=True)
    scratch_work.mkdir(parents=True, exist_ok=True)

    reports_dir = project_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    log_file = reports_dir / f"nextflow_{entrypoint}.log"

    # Pre-generate unrolled DAG pipeline snapshot
    total_ns = 60.0
    chunk_ns = 10.0
    if configs:
        prod_cfg = configs.get("md_simulation", {}).get("production", {})
        total_ns = float(prod_cfg.get("total_time_ns", 60.0))
        chunk_ns = float(prod_cfg.get("chunk_size_ns", 10.0))

    for i, arg in enumerate(extra_args):
        if arg == "--md_total_ns" and i + 1 < len(extra_args):
            total_ns = float(extra_args[i + 1])
        elif arg == "--md_chunk_ns" and i + 1 < len(extra_args):
            chunk_ns = float(extra_args[i + 1])

    total_c = max(1, int(round(total_ns / chunk_ns)))

    pipeline_code_dir = scratch_base / ".pipeline_code"
    if pipeline_code_dir.exists():
        shutil.rmtree(pipeline_code_dir)
    pipeline_code_dir.mkdir(parents=True)

    for d in ["workflows", "subworkflows", "modules", "assets", "conf", "bin"]:
        src = FRAMEWORK_DIR / d
        if src.exists():
            shutil.copytree(src, pipeline_code_dir / d)
    for f in ["main.nf", "nextflow.config"]:
        src = FRAMEWORK_DIR / f
        if src.exists():
            shutil.copy(src, pipeline_code_dir / f)

    unrolled_code = build_unrolled_md_dag(total_c)
    (pipeline_code_dir / "subworkflows" / "run_chunked_md.nf").write_text(unrolled_code)

    main_nf = pipeline_code_dir / "main.nf"

    for script_file in (pipeline_code_dir / "bin").glob("*.py"):
        try:
            script_file.chmod(0o755)
        except Exception:
            pass

    env = os.environ.copy()
    env["NXF_HOME"] = str(scratch_nxf_home.resolve())
    env["NXF_WORK"] = str(scratch_work.resolve())

    cmd = [
        "nextflow",
        "-log",
        str(log_file.resolve()),
        "run",
        str(main_nf.resolve()),
        "-work-dir",
        str(scratch_work.resolve()),
        "-profile",
        profile,
        "--entry",
        entrypoint,
    ]

    if resume:
        cmd.append("-resume")

    cmd.extend(extra_args)

    console.rule(
        f"[bold cyan]Nextflow Orchestrator: {entrypoint.upper()} ({profile.upper()})[/bold cyan]"
    )
    console.print(f"[dim]Project Scratch: {scratch_base}[/dim]")
    console.print(f"[dim]Deliverables:    {project_dir.resolve()}[/dim]")
    console.print(f"[bold blue]Command:[/bold blue] {' '.join(cmd)}\n")

    # Run Nextflow INSIDE the scratch directory.
    # This prevents .nextflow/ caches and .nextflow.log from polluting the user's clean project folder.
    try:
        proc = subprocess.Popen(cmd, cwd=scratch_base, env=env)
        ret = proc.wait()

        if ret == 0:
            console.print(
                f"\n[bold green]✓ Subworkflow '{entrypoint}' finished successfully![/bold green]"
            )
            if auto_clean:
                clean_shared_scratch(str(scratch_base), verbose=True)
            return 0
        else:
            console.print(
                f"\n[bold red]✗ Nextflow failed with exit code {ret}[/bold red]"
            )
            console.print(
                f"[yellow]Scratch directory {scratch_base} preserved for '-resume'.[/yellow]"
            )
            return ret

    except KeyboardInterrupt:
        console.print(
            "\n[bold yellow]Execution canceled by user. Scratch preserved for resume.[/bold yellow]"
        )
        return 130
    except Exception as e:
        console.print(f"\n[bold red]Failed to launch Nextflow: {e}[/bold red]")
        return 1
