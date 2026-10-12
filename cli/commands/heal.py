"""
File: cli/commands/heal.py
Description: Resolves broken crystal structures using physical AI inpainting.
"""

import json
import sys
from pathlib import Path
from rich.console import Console
from rich.prompt import Prompt
from cli.core.db import init_db, Target, ReceptorFrame
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()

def cmd_heal(args):
    cwd = Path.cwd()
    session = init_db()

    query = session.query(ReceptorFrame).filter_by(state="broken")
    if getattr(args, "target", None):
        t = session.query(Target).filter_by(name=args.target).first()
        if not t:
            console.print(f"[bold red]FATAL: Target '{args.target}' not found.[/bold red]")
            sys.exit(1)
        query = query.filter_by(target_id=t.id)

    broken_frames = query.all()
    if not broken_frames:
        console.print("[bold green]✓ No broken frames found. Everything is whole and ready to prep![/bold green]")
        session.close()
        return

    # Check for missing FASTA sequences and prompt interactively
    targets_needing_fasta = set(f.target for f in broken_frames if not f.target.reference_sequence)
    for t in targets_needing_fasta:
        console.print(f"\n[bold yellow]Target '{t.name}' has broken frames but NO reference FASTA sequence.[/bold yellow]")
        console.print("Paste the complete FASTA sequence below. Press [bold]Enter[/bold], then [bold]Ctrl+D[/bold] to submit:")
        try:
            seq_input = sys.stdin.read()
        except EOFError:
            pass
        
        clean_seq = "".join([line.strip() for line in seq_input.splitlines() if not line.startswith(">")])
        if not clean_seq:
            console.print("[bold red]FATAL: Empty sequence provided. Cannot heal without reference.[/bold red]")
            sys.exit(1)
            
        t.reference_sequence = clean_seq.upper()
        session.commit()
        console.print(f"[green]✓ Sequence saved for Target '{t.name}'[/green]")

    heal_plan = []
    for f in broken_frames:
        heal_plan.append({
            "target_name": f.target.name,
            "frame_id": f.name,
            "path_raw": str(Path(f.path).resolve()),
            "fasta_seq": f.target.reference_sequence
        })

    plan_file = cwd / "heal_plan.json"
    with open(plan_file, "w", encoding="utf-8") as pf:
        json.dump(heal_plan, pf, indent=2)

    console.print(f"[bold cyan]Dispatching {len(heal_plan)} frame(s) to AI surgical healing module...[/bold cyan]")

    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="heal")
    nf_args.extend(["--heal_plan", str(plan_file.resolve())])

    if getattr(args, "urgent", False):
        nf_args.extend(["--queue", "urgent"])

    profile = "local" if getattr(args, "local", False) else "slurm"
    resume = not getattr(args, "no_resume", False)

    ret_code = execute_nextflow("heal", nf_args, cwd, profile=profile, resume=resume, configs=configs)

    if ret_code != 0:
        console.print("[bold red]Nextflow healing failed.[/bold red]")
        sys.exit(ret_code)

    receipt = cwd / "heal_receipt.tsv"
    if not receipt.exists() or receipt.stat().st_size == 0:
        console.print("[bold red]FATAL: Execution receipt not found. Pipeline failed silently.[/bold red]")
        sys.exit(1)

    healed_count = 0
    with open(receipt, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip(): continue
            frame_name, healed_path = line.strip().split("\t")
            h_path = Path(healed_path).resolve()
            
            if h_path.exists() and h_path.stat().st_size > 0:
                frame = session.query(ReceptorFrame).filter_by(name=frame_name).first()
                if frame:
                    frame.path = str(h_path)
                    frame.state = "raw"
                    frame.source_type = "crystal_healed"
                    healed_count += 1

    session.commit()
    session.close()

    console.print(f"\n[bold green]✓ Successfully healed {healed_count} frame(s)! They are now in 'raw' state.[/bold green]")
    console.print("[dim]Use 'sbdd prep' to process them for docking or MD.[/dim]")
