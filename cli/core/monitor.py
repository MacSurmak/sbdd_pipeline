import socket
import json
import re
from rich.console import Console
from rich.progress import (
    Progress,
    BarColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

console = Console()


def run_monitor(port: int = 8089, default_chunk_steps: int = 5000000):
    """Listens for UDP telemetry packets and renders real-time progress bars."""
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.bind(("0.0.0.0", port))
    sock.settimeout(0.5)

    console.rule("[bold cyan]SBDD Cluster MD & Docking Live Monitor[/bold cyan]")
    console.print(
        f"[dim]Listening for UDP telemetry on port {port}... Press Ctrl+C to exit.[/dim]\n"
    )

    progress = Progress(
        TextColumn("[bold cyan]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TextColumn("({task.completed}/{task.total})"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
    )

    task_bars = {}

    try:
        with progress:
            while True:
                try:
                    data, _ = sock.recvfrom(2048)
                    msg = json.loads(data.decode("utf-8", errors="ignore"))
                    task_name = msg.get("task", "unknown")
                    log_line = msg.get("log", "")

                    # 1. GROMACS Step Telemetry
                    # 1. GROMACS Step Telemetry (Matches 'Step 500' and 'Step= 50')
                    match_gmx = re.search(r"(?i)step[=\s]+(\d+)", log_line)
                    if match_gmx:
                        current_step = int(match_gmx.group(1))
                        if task_name not in task_bars:
                            # Strip Nextflow pipeline paths for clean display
                            clean_name = task_name.split(":")[-1].strip()
                            clean_name = clean_name.replace(
                                "RUN_MD_CHUNK ", ""
                            ).replace("GMX_EQUILIBRATION ", "EQ: ")
                            is_eq = "EQUI" in task_name.upper()
                            total = 50000 if is_eq else default_chunk_steps
                            task_bars[task_name] = progress.add_task(
                                f"[green]{clean_name}", total=total
                            )
                        progress.update(task_bars[task_name], completed=current_step)
                        continue

                    # 2. GNINA Percentage Telemetry
                    match_gnina = re.search(r"\[\s*(\d+)%\]", log_line)
                    if match_gnina:
                        pct = int(match_gnina.group(1))
                        if task_name not in task_bars:
                            clean_name = task_name.replace(
                                "GNINA_DOCK ", "DOCK: "
                            ).replace("GNINA_RESCORE ", "RESC: ")
                            task_bars[task_name] = progress.add_task(
                                f"[magenta]{clean_name}", total=100
                            )
                        progress.update(task_bars[task_name], completed=pct)

                except socket.timeout:
                    continue
                except Exception:
                    pass
    except KeyboardInterrupt:
        console.print("\n[yellow]Monitoring stopped.[/yellow]")
    finally:
        sock.close()
