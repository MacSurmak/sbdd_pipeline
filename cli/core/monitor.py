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


def run_monitor(port: int = 8088):
    """
    Renders high-frequency (10Hz) real-time cluster progress bars by listening to the local
    UDP fanout stream on 127.0.0.1:8088 relayed by sbdd-telemetry.
    Zero-delay updates with smooth GPU power and utilization metrics.
    """
    console.rule("[bold cyan]SBDD Cluster MD Live Monitor (GROMACS)[/bold cyan]")
    console.print(
        f"[dim]Streaming live 10Hz telemetry from 127.0.0.1:{port}... Press Ctrl+C to exit.[/dim]\n"
    )

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    if hasattr(socket, "SO_REUSEPORT"):
        try:
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEPORT, 1)
        except Exception:
            pass
    sock.bind(("127.0.0.1", port))
    sock.settimeout(0.2)

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
    last_seen = {}

    try:
        with progress:
            while True:
                now = time.time()

                # Remove stalled or halted tasks after 3.5s of silence
                for t_id in list(task_bars.keys()):
                    if now - last_seen.get(t_id, now) > 3.5:
                        progress.remove_task(task_bars[t_id])
                        del task_bars[t_id]
                        last_seen.pop(t_id, None)

                try:
                    data, _ = sock.recvfrom(65535)
                    payload = json.loads(data.decode("utf-8", errors="ignore"))
                    t_id = payload.get("task", "unknown")
                    if t_id == "unknown":
                        continue

                    if payload.get("event") == "exit":
                        if t_id in task_bars:
                            progress.remove_task(task_bars[t_id])
                            del task_bars[t_id]
                            last_seen.pop(t_id, None)
                        continue

                    last_seen[t_id] = now

                    step = payload.get("display_step", payload.get("step")) or 0
                    total = payload.get("display_total", payload.get("total")) or 1
                    dt_ps = payload.get("dt_ps") or 0.002

                    step_info = payload.get("step_info")
                    if not step_info:
                        t_str = f" [dim]({format_time_ns(step * dt_ps)})[/dim]"
                        step_info = f"[cyan]{format_steps(step)}/{format_steps(total)}[/cyan]{t_str}"

                    speed_str = payload.get("speed", "-- ns/d")
                    eta_str = payload.get("eta", "--:--")

                    hw_gpu = payload.get("hardware", {}).get("gpu", {})
                    gpu_name = hw_gpu.get("model") or "GPU"
                    if "RTX" in gpu_name:
                        gpu_short = "RTX " + gpu_name.split("RTX")[-1].strip()
                    else:
                        gpu_short = gpu_name

                    gpu_util = hw_gpu.get("utilization_gpu_pct")
                    gpu_str = f"{gpu_short} ({gpu_util}%)" if gpu_util is not None else gpu_short

                    p_watts = hw_gpu.get("power_watts")
                    power_str = f"{p_watts:>5.1f} W" if p_watts is not None else "-- W"

                    clean_name = t_id.strip()
                    m_chunk_wrap = re.match(
                        r"^RUN_MD_CHUNK_\d+_(.+?)\s+(\(#\d+\))\s*\([^)]*\)$", clean_name
                    )
                    if m_chunk_wrap:
                        clean_name = f"MD: {m_chunk_wrap.group(1)} {m_chunk_wrap.group(2)}"
                    elif ":" in clean_name and not clean_name.startswith(("EQ:", "MD:", "SWISH:")):
                        clean_name = clean_name.split(":")[-1].strip()

                    node_str = payload.get("node") or "node"

                    if t_id not in task_bars:
                        task_bars[t_id] = progress.add_task(
                            f"[green]{escape(clean_name)}",
                            total=total,
                            completed=step,
                            node=escape(node_str),
                            gpu=escape(gpu_str),
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
                            node=escape(node_str),
                            gpu=escape(gpu_str),
                            power=power_str,
                            speed=speed_str,
                            step_info=step_info,
                            eta=eta_str,
                        )
                except socket.timeout:
                    continue
                except Exception:
                    pass

    except KeyboardInterrupt:
        console.print("\n[yellow]Monitoring stopped.[/yellow]")
    finally:
        sock.close()
