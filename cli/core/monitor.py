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


from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn, TimeRemainingColumn, MofNCompleteColumn

console = Console()

def run_monitor(port: int = 8089):
    """Listens for UDP telemetry packets and renders real-time progress bars with dynamic nsteps."""
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.bind(("0.0.0.0", port))
    sock.settimeout(0.5)

    console.rule("[bold cyan]SBDD Cluster MD & Docking Live Monitor[/bold cyan]")
    console.print(f"[dim]Listening for UDP telemetry on port {port}... Press Ctrl+C to exit.[/dim]\n")

    progress = Progress(
        TextColumn("[bold cyan]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console
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

                    clean_name = task_name.split(":")[-1].strip()
                    clean_name = clean_name.replace("RUN_MD_CHUNK ", "").replace("GMX_EQUILIBRATION ", "EQ: ")

                    # 1. Dynamically parse target total steps from GROMACS startup header
                    match_total = re.search(r'(\d+)\s+steps,\s+[\d.]+\s+ps', log_line) or \
                                  re.search(r'(?i)number of steps\s*=\s*(\d+)', log_line) or \
                                  re.search(r'(?i)overriding nsteps.*:\s*(\d+)\s+steps', log_line)
                    if match_total:
                        total_steps = int(match_total.group(1))
                        if task_name not in task_bars:
                            task_bars[task_name] = progress.add_task(f"[green]{clean_name}", total=total_steps)
                        else:
                            progress.update(task_bars[task_name], total=total_steps)
                        continue

                    # 2. GROMACS Step Progress (Matches 'step 11280500' and 'Step= 50')
                    match_step = re.search(r'(?i)step[=\s]+(\d+)', log_line)
                    if match_step:
                        current_step = int(match_step.group(1))
                        if task_name not in task_bars:
                            # Total is unknown until header arrives; renders indeterminate pulse/counter cleanly
                            task_bars[task_name] = progress.add_task(f"[green]{clean_name}", total=None)

                        # Auto-expand total if run exceeded previous estimate
                        task_obj = progress.tasks[task_bars[task_name]]
                        if task_obj.total is not None and current_step > task_obj.total:
                            progress.update(task_bars[task_name], total=None)

                        progress.update(task_bars[task_name], completed=current_step)
                        continue

                    # 3. GNINA Percentage Telemetry
                    match_gnina = re.search(r'\[\s*(\d+)%\]', log_line)
                    if match_gnina:
                        pct = int(match_gnina.group(1))
                        if task_name not in task_bars:
                            clean_name = task_name.replace("GNINA_DOCK ", "DOCK: ").replace("GNINA_RESCORE ", "RESC: ")
                            task_bars[task_name] = progress.add_task(f"[magenta]{clean_name}", total=100)
                        progress.update(task_bars[task_name], completed=pct)                            ).replace("GNINA_RESCORE ", "RESC: ")
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
