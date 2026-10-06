from datetime import datetime, timedelta
import json
import math
import os
import re
import socket
import subprocess
import time
from rich.console import Console
from rich.markup import escape
from rich.progress import (
    BarColumn,
    Progress,
    TextColumn,
    TimeElapsedColumn,
)

console = Console()


def format_steps(n: int | float | None) -> str:
    """Format step counts into dynamic readable metric units (units, K, M, G)."""
    if n is None or not math.isfinite(float(n)):
        return "?"
    val = float(n)
    abs_val = abs(val)
    if abs_val >= 1e9:
        return f"{val / 1e9:.2f}G"
    elif abs_val >= 1e6:
        return f"{val / 1e6:.2f}M"
    elif abs_val >= 1e3:
        if abs_val >= 100e3:
            return f"{val / 1e3:.0f}K"
        return f"{val / 1e3:.1f}K"
    else:
        return str(int(round(val)))


def format_time_ns(ps: float | None) -> str:
    """Format picoseconds into physical simulation time (ps, ns, µs)."""
    if ps is None or not math.isfinite(float(ps)):
        return ""
    ns = ps / 1000.0
    abs_ns = abs(ns)
    if abs_ns >= 1000.0:
        return f"{ns / 1000.0:.2f} µs"
    elif abs_ns >= 1.0:
        return f"{ns:.1f} ns"
    else:
        return f"{ps:.1f} ps"


def format_eta(seconds_left: float | None) -> str:
    """Format remaining seconds into human-readable duration and projected calendar completion date/time."""
    if seconds_left is None or not math.isfinite(seconds_left) or seconds_left < 0:
        return "--:--"

    now = datetime.now()
    finish_dt = now + timedelta(seconds=seconds_left)

    if seconds_left >= 86400:
        days = int(seconds_left // 86400)
        hours = int((seconds_left % 86400) // 3600)
        dur_str = f"{days}d {hours:02d}h"
    elif seconds_left >= 3600:
        hours = int(seconds_left // 3600)
        mins = int((seconds_left % 3600) // 60)
        dur_str = f"{hours}h {mins:02d}m"
    elif seconds_left >= 60:
        mins = int(seconds_left // 60)
        secs = int(seconds_left % 60)
        dur_str = f"{mins}m {secs:02d}s"
    else:
        dur_str = f"{int(seconds_left)}s"

    if finish_dt.year != now.year:
        date_str = finish_dt.strftime("%d.%m.%Y %H:%M")
    elif (finish_dt.date() - now.date()).days > 0:
        date_str = finish_dt.strftime("%d.%m %H:%M")
    else:
        date_str = finish_dt.strftime("%H:%M")

    return f"{dur_str} (ETA: {date_str})"


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
        TextColumn("{task.fields[step_info]}"),
        TextColumn("[bold green]{task.fields[speed]}[/bold green]"),
        TimeElapsedColumn(),
        TextColumn("[yellow]{task.fields[eta]}[/yellow]"),
        console=console,
    )

    task_bars = {}
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
                        # When job has a Slurm ID and is no longer present in squeue
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

                    # Handle explicit exit signal from worker telemetry emitter
                    if event == "exit":
                        if task_name in task_bars:
                            progress.remove_task(task_bars[task_name])
                            del task_bars[task_name]
                            last_seen.pop(task_name, None)
                        continue

                    last_seen[task_name] = now

                    # Clean task description and normalize container wrapper task tags
                    clean_name = task_name.strip()
                    m_chunk_wrap = re.match(
                        r"^RUN_MD_CHUNK_\d+_(.+?)\s+(\(#\d+\))\s*\([^)]*\)$", clean_name
                    )
                    if m_chunk_wrap:
                        clean_name = (
                            f"MD: {m_chunk_wrap.group(1)} {m_chunk_wrap.group(2)}"
                        )
                    elif ":" in clean_name and not clean_name.startswith(
                        ("EQ:", "MD:", "SWISH:")
                    ):
                        clean_name = clean_name.split(":")[-1].strip()

                    node_name = msg.get("node", "local")

                    # Extract authoritative pre-computed metrics emitted by worker
                    display_step = msg.get("display_step", msg.get("step"))
                    display_total = msg.get("display_total", msg.get("total"))
                    step_info_str = msg.get("step_info")
                    speed_display = msg.get("speed", "-- ns/day")
                    eta_display = msg.get("eta", "--:--")

                    # Fallback formatting if raw telemetry packet arrives without pre-rendered strings
                    if step_info_str is None and display_step is not None:
                        dt_ps = msg.get("dt_ps")
                        t_str = f" [dim]({format_time_ns(display_step * dt_ps)})[/dim]" if dt_ps else ""
                        if display_total is not None:
                            step_info_str = f"[cyan]{format_steps(display_step)}/{format_steps(display_total)}[/cyan]{t_str}"
                        else:
                            step_info_str = f"[cyan]{format_steps(display_step)}[/cyan]{t_str}"

                    if display_step is not None:
                        if task_name not in task_bars:
                            task_bars[task_name] = progress.add_task(
                                f"[green]{escape(clean_name)}",
                                total=display_total,
                                completed=display_step,
                                node=escape(node_name),
                                speed=speed_display,
                                step_info=step_info_str or "",
                                eta=eta_display,
                            )
                        else:
                            task_obj = progress._tasks[task_bars[task_name]]
                            if (
                                display_total is not None
                                and task_obj.total != display_total
                            ):
                                progress.update(
                                    task_bars[task_name], total=display_total
                                )

                            progress.update(
                                task_bars[task_name],
                                completed=display_step,
                                node=escape(node_name),
                                speed=speed_display,
                                step_info=step_info_str or "",
                                eta=eta_display,
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
