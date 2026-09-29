import socket
import json
import math
import re
import time
import subprocess
from rich.console import Console
from rich.markup import escape
from rich.progress import (
    Progress,
    BarColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
    MofNCompleteColumn,
)

console = Console()


def run_monitor(port: int = 8089):
    """Listens for UDP telemetry packets from GROMACS and renders real-time progress bars with dynamic lifecycle."""
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.bind(("0.0.0.0", port))
    sock.settimeout(0.5)

    console.rule("[bold cyan]SBDD Cluster MD Live Monitor (GROMACS)[/bold cyan]")
    console.print(
        f"[dim]Listening for GROMACS UDP telemetry on port {port}... Press Ctrl+C to exit.[/dim]\n"
    )

    progress = Progress(
        TextColumn("[bold cyan]{task.description}"),
        TextColumn("[bold yellow]({task.fields[node]})[/bold yellow]"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        MofNCompleteColumn(),
        TextColumn("[bold green]{task.fields[speed]}[/bold green]"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
    )

    task_bars = {}
    task_stats = (
        {}
    )  # task_name -> {"last_step": int, "last_time": float, "ema_rate": float, "dt_ps": float, "speed_str": str}
    last_seen = {}
    last_reap = time.time()

    try:
        with progress:
            while True:
                now = time.time()

                # Periodic reaper: sync with squeue to purge jobs terminated via SIGKILL, scancel or completion
                if now - last_reap > 2.0:
                    last_reap = now
                    active_job_ids = None
                    try:
                        sq_out = subprocess.check_output(
                            ["squeue", "-h", "-o", "%i"],
                            text=True,
                            stderr=subprocess.DEVNULL,
                        )
                        active_job_ids = set(sq_out.split())
                    except Exception:
                        pass

                    for t_name in list(task_bars.keys()):
                        m_job = re.search(r"#(\d+)", t_name)
                        # When job has a Slurm ID and is no longer present in squeue (even if squeue is completely empty)
                        if (
                            m_job
                            and active_job_ids is not None
                            and m_job.group(1) not in active_job_ids
                        ):
                            progress.remove_task(task_bars[t_name])
                            del task_bars[t_name]
                            last_seen.pop(t_name, None)
                            continue

                        # If task went completely silent for over 25s
                        if now - last_seen.get(t_name, now) > 25.0:
                            progress.remove_task(task_bars[t_name])
                            del task_bars[t_name]
                            last_seen.pop(t_name, None)

                try:
                    data, _ = sock.recvfrom(2048)
                    msg = json.loads(data.decode("utf-8", errors="ignore"))
                    task_name = msg.get("task", "unknown")
                    event = msg.get("event")

                    # Handle explicit exit signal from wrapper trap
                    if event == "exit":
                        if task_name in task_bars:
                            progress.remove_task(task_bars[task_name])
                            del task_bars[task_name]
                            last_seen.pop(task_name, None)
                        continue

                    last_seen[task_name] = now
                    log_line = msg.get("log", "").replace("\r", " ").strip()

                    # Make description compact so it never wraps on narrow terminal windows
                    clean_name = task_name.split(":")[-1].strip()
                    clean_name = clean_name.replace("GMX_EQUILIBRATION_", "EQ: ")
                    clean_name = clean_name.replace("GMX_EQUILIBRATION", "EQ")
                    clean_name = clean_name.replace("RUN_MD_CHUNK_", "MD: ")
                    clean_name = clean_name.replace("RUN_MD_CHUNK", "MD")

                    node_name = msg.get("node", "local")
                    dt_ps = msg.get("dt_ps")

                    # Extract total from packet metadata or log line
                    total_steps = msg.get("total")
                    match_total = re.search(r"(\d+)\s+steps,\s+([\d.]+)\s+ps", log_line)
                    if match_total:
                        total_steps = int(match_total.group(1))
                        total_ps = float(match_total.group(2))
                        if total_steps > 0:
                            dt_ps = total_ps / total_steps
                    else:
                        m_steps_only = re.search(
                            r"(?i)number of steps\s*=\s*(\d+)", log_line
                        ) or re.search(
                            r"(?i)overriding nsteps.*:\s*(\d+)\s+steps", log_line
                        )
                        if m_steps_only:
                            total_steps = int(m_steps_only.group(1))

                    # GROMACS Step Progress (Matches 'step 11280500' and 'Step= 50')
                    match_step = re.search(r"(?i)step[=\s]+(\d+)", log_line)
                    if match_step:
                        current_step = int(match_step.group(1))

                        if task_name not in task_stats:
                            task_stats[task_name] = {
                                "last_step": current_step,
                                "last_time": now,
                                "ema_rate": None,
                                "dt_ps": dt_ps,
                                "speed_str": "-- ns/day" if dt_ps else "-- step/s",
                            }
                        else:
                            stats = task_stats[task_name]
                            if dt_ps and not stats["dt_ps"]:
                                stats["dt_ps"] = dt_ps

                            # Reset on stage transition or step rollback
                            if current_step < stats["last_step"]:
                                stats["last_step"] = current_step
                                stats["last_time"] = now
                                stats["ema_rate"] = None

                            dt_wall = now - stats["last_time"]
                            if dt_wall >= 1.0:
                                d_step = current_step - stats["last_step"]
                                if d_step >= 0:
                                    inst_rate = d_step / dt_wall
                                    if stats["ema_rate"] is None:
                                        stats["ema_rate"] = inst_rate
                                    else:
                                        # Exponential moving average with time constant tau = 8.0s
                                        alpha = 1.0 - math.exp(-dt_wall / 8.0)
                                        stats["ema_rate"] = (
                                            alpha * inst_rate
                                            + (1.0 - alpha) * stats["ema_rate"]
                                        )

                                    smoothed_rate = stats["ema_rate"]
                                    if stats["dt_ps"]:
                                        ns_day = smoothed_rate * stats["dt_ps"] * 86.4
                                        stats["speed_str"] = f"{ns_day:>5.1f} ns/d"
                                    else:
                                        stats["speed_str"] = (
                                            f"{smoothed_rate:>5.0f} st/s"
                                        )

                                stats["last_step"] = current_step
                                stats["last_time"] = now

                        speed_display = task_stats[task_name]["speed_str"]

                        if task_name not in task_bars:
                            task_bars[task_name] = progress.add_task(
                                f"[green]{escape(clean_name)}",
                                total=total_steps,
                                node=escape(node_name),
                                speed=speed_display,
                            )
                        else:
                            task_obj = progress._tasks[task_bars[task_name]]
                            if (
                                total_steps is not None
                                and task_obj.total != total_steps
                            ):
                                progress.update(task_bars[task_name], total=total_steps)

                            progress.update(
                                task_bars[task_name],
                                completed=current_step,
                                node=escape(node_name),
                                speed=speed_display,
                            )
                        continue

                except socket.timeout:
                    continue
                except Exception:
                    pass
    except KeyboardInterrupt:
        console.print("\n[yellow]Monitoring stopped.[/yellow]")
    finally:
        sock.close()
