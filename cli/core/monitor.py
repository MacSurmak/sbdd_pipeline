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


import sqlite3
from pathlib import Path

TELEMETRY_DB = Path("/var/lib/sbdd/cluster_telemetry.db")


def run_monitor(port: int = 8089):
    """
    Renders real-time cluster progress bars by streaming active simulations directly
    from the local telemetry WAL SQLite database maintained by sbdd-telemetry.
    Eliminates socket address collisions and displays live GPU wattage and hardware load.
    """
    console.rule("[bold cyan]SBDD Cluster MD Live Monitor (GROMACS)[/bold cyan]")
    console.print(
        f"[dim]Streaming live telemetry from {TELEMETRY_DB}... Press Ctrl+C to exit.[/dim]\n"
    )

    progress = Progress(
        TextColumn("[bold cyan]{task.description}"),
        TextColumn("[bold yellow]({task.fields[node]})[/bold yellow]"),
        TextColumn("[magenta]{task.fields[gpu]}[/magenta]"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TextColumn("{task.fields[step_info]}"),
        TextColumn("[bold green]{task.fields[speed]}[/bold green]"),
        TextColumn("[yellow]{task.fields[power]}[/yellow]"),
        TextColumn("[dim]{task.fields[eta]}[/dim]"),
        console=console,
    )

    task_bars = {}

    try:
        with progress:
            while True:
                if not TELEMETRY_DB.exists():
                    time.sleep(1.0)
                    continue

                active_tasks = {}
                try:
                    conn = sqlite3.connect(
                        f"file:{TELEMETRY_DB}?mode=ro", uri=True, timeout=1.0
                    )
                    conn.row_factory = sqlite3.Row
                    cursor = conn.cursor()
                    cursor.execute(
                        """
                        SELECT
                            task_id, node_name, stage, natoms, timestep_fs,
                            total_steps, last_step, speed_cumulative_ns_day,
                            speed_instant_ns_day, eta_string, gpu_model,
                            gpu_power_avg_w, gpu_vram_peak_mb, status,
                            strftime('%s', 'now') - strftime('%s', updated_at) AS age_sec
                        FROM simulation_runs
                        WHERE status = 'running'
                        ORDER BY updated_at DESC;
                        """
                    )
                    rows = cursor.fetchall()
                    conn.close()

                    for r in rows:
                        active_tasks[r["task_id"]] = dict(r)
                except Exception:
                    pass

                # Purge bars for tasks that are no longer actively running
                for t_name in list(task_bars.keys()):
                    if t_name not in active_tasks:
                        progress.remove_task(task_bars[t_name])
                        del task_bars[t_name]

                # Update or add progress bars
                for t_id, data in active_tasks.items():
                    step = data.get("last_step") or 0
                    total = data.get("total_steps") or 1
                    dt_ps = (data.get("timestep_fs") or 2.0) / 1000.0

                    t_str = f" [dim]({format_time_ns(step * dt_ps)})[/dim]"
                    step_info = f"[cyan]{format_steps(step)}/{format_steps(total)}[/cyan]{t_str}"

                    speed_val = data.get("speed_cumulative_ns_day") or data.get("speed_instant_ns_day")
                    speed_str = f"{speed_val:>5.1f} ns/d" if speed_val else "-- ns/d"

                    gpu_short = data.get("gpu_model") or "GPU"
                    if "RTX" in gpu_short:
                        gpu_short = "RTX " + gpu_short.split("RTX")[-1].strip()

                    power_val = data.get("gpu_power_avg_w")
                    power_str = f"{power_val:>5.1f} W" if power_val else "-- W"

                    eta_str = data.get("eta_string") or "--:--"

                    # Normalize title display
                    clean_name = t_id.strip()
                    m_chunk_wrap = re.match(
                        r"^RUN_MD_CHUNK_\d+_(.+?)\s+(\(#\d+\))\s*\([^)]*\)$", clean_name
                    )
                    if m_chunk_wrap:
                        clean_name = f"MD: {m_chunk_wrap.group(1)} {m_chunk_wrap.group(2)}"
                    elif ":" in clean_name and not clean_name.startswith(("EQ:", "MD:", "SWISH:")):
                        clean_name = clean_name.split(":")[-1].strip()

                    if t_id not in task_bars:
                        task_bars[t_id] = progress.add_task(
                            f"[green]{escape(clean_name)}",
                            total=total,
                            completed=step,
                            node=escape(data.get("node_name") or "node"),
                            gpu=escape(gpu_short),
                            power=power_str,
                            speed=speed_str,
                            step_info=step_info,
                            eta=eta_str,
                        )
                    else:
                        task_obj = progress._tasks[task_bars[t_id]]
                        if total is not None and task_obj.total != total:
                            progress.update(task_bars[t_id], total=total)

                        progress.update(
                            task_bars[t_id],
                            completed=step,
                            node=escape(data.get("node_name") or "node"),
                            gpu=escape(gpu_short),
                            power=power_str,
                            speed=speed_str,
                            step_info=step_info,
                            eta=eta_str,
                        )

                time.sleep(0.5)

    except KeyboardInterrupt:
        console.print("\n[yellow]Monitoring stopped.[/yellow]")
