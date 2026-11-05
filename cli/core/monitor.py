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
    Renders high-frequency (20 FPS) real-time cluster progress bars by listening to the local
    UDP fanout stream on 127.0.0.1:8088 relayed by sbdd-telemetry.
    Monitors all atomic cluster tasks: GROMACS MD, Glide, GNINA, Prime, MM-PBSA, Packmol, Declash, and ORCA.
    """
    console.rule(
        "[bold cyan]SBDD Cluster Live Monitor (MD / Docking / QM / Analysis)[/bold cyan]"
    )
    console.print(
        f"[dim]Live telemetry from 127.0.0.1:{port}... Press Ctrl+C to exit.[/dim]\n"
    )

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    if hasattr(socket, "SO_REUSEPORT"):
        try:
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEPORT, 1)
        except Exception:
            pass
    sock.bind(("127.0.0.1", port))
    sock.settimeout(0.05)

    progress = Progress(
        TextColumn("[bold cyan]{task.description}"),
        TextColumn("[bold yellow]{task.fields[node_display]}[/bold yellow]"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TextColumn("{task.fields[step_info]}"),
        TextColumn("[bold green]{task.fields[speed]}[/bold green]"),
        TextColumn("[magenta]{task.fields[gpu]}[/magenta]"),
        TextColumn("[yellow]{task.fields[power]}[/yellow]"),
        TextColumn("[dim]{task.fields[eta]}[/dim]"),
        console=console,
    )

    task_bars = {}
    last_seen = {}

    task_real_step = {}
    task_packet_time = {}
    task_rate_steps_s = {}
    task_totals = {}
    task_dts = {}
    task_raw_step_info = {}
    task_nodes = {}
    task_gpus = {}
    task_powers = {}
    task_speeds = {}
    task_etas = {}
    task_names = {}

    task_parents = {}
    task_children = {}
    task_is_subtask = {}
    task_subtask_label = {}

    cpu_only_tools = {
        "glide",
        "prime",
        "prepwizard",
        "orca",
        "mmpbsa",
        "packmol",
        "declash",
    }

    try:
        with progress:
            while True:
                now = time.time()

                while True:
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
                            p_id = task_parents.pop(t_id, None)
                            if p_id and p_id in task_children:
                                task_children[p_id].discard(t_id)
                            for child_id in list(task_children.pop(t_id, set())):
                                if child_id in task_bars:
                                    progress.remove_task(task_bars[child_id])
                                    del task_bars[child_id]
                                task_parents.pop(child_id, None)
                                last_seen.pop(child_id, None)
                                task_real_step.pop(child_id, None)

                            last_seen.pop(t_id, None)
                            task_real_step.pop(t_id, None)
                            continue

                        last_seen[t_id] = now
                        tool_raw = (payload.get("tool") or "mdrun").lower()
                        parent_id = payload.get("parent_task")
                        is_sub = bool(payload.get("is_subtask") or parent_id)
                        task_is_subtask[t_id] = is_sub

                        if is_sub and parent_id:
                            task_parents[t_id] = parent_id
                            task_children.setdefault(parent_id, set()).add(t_id)
                            task_subtask_label[t_id] = payload.get("subtask_label") or t_id

                        step = payload.get("display_step", payload.get("step")) or 0
                        total = payload.get("display_total", payload.get("total")) or 1
                        dt_ps = payload.get("dt_ps")

                        task_real_step[t_id] = step
                        task_packet_time[t_id] = now
                        task_totals[t_id] = total
                        task_dts[t_id] = dt_ps
                        task_raw_step_info[t_id] = payload.get("step_info")

                        perf = payload.get("performance", {})
                        rep_speed = perf.get("instant_ns_day") or perf.get(
                            "cumulative_ns_day"
                        )
                        if rep_speed and dt_ps and dt_ps > 0:
                            task_rate_steps_s[t_id] = rep_speed / (dt_ps * 86.4)
                        else:
                            task_rate_steps_s[t_id] = 0.0

                        hw = payload.get("hardware", {})
                        hw_gpu = hw.get("gpu", {})
                        hw_cpu = hw.get("cpu", {})
                        hw_mem = hw.get("memory", {})

                        if tool_raw in cpu_only_tools:
                            cores = hw_cpu.get("allocated_cores") or 1
                            task_gpus[t_id] = f"CPU ({cores}c)"
                            ram_used = hw_mem.get("ram_used_mb")
                            task_powers[t_id] = (
                                f"{ram_used / 1024.0:>4.1f} GB"
                                if ram_used
                                else "-- W"
                            )
                        else:
                            gpu_name = hw_gpu.get("model") or "GPU"
                            if "RTX" in gpu_name:
                                gpu_short = "RTX " + gpu_name.split("RTX")[-1].strip()
                            elif (
                                "Tesla" in gpu_name
                                or "A100" in gpu_name
                                or "H100" in gpu_name
                            ):
                                gpu_short = gpu_name.split()[-1]
                            else:
                                gpu_short = gpu_name

                            gpu_util = hw_gpu.get("utilization_gpu_pct")
                            task_gpus[t_id] = (
                                f"{gpu_short} ({gpu_util}%)"
                                if gpu_util is not None
                                else gpu_short
                            )

                            p_watts = hw_gpu.get("power_watts")
                            task_powers[t_id] = (
                                f"{p_watts:>5.1f} W" if p_watts is not None else "-- W"
                            )

                        multisim = payload.get("simulation", {}).get("multisim", {})
                        n_reps = (
                            int(multisim.get("replicas") or 1)
                            if multisim.get("enabled")
                            else 1
                        )

                        if n_reps > 1 and rep_speed:
                            agg_speed = perf.get("aggregate_instant_ns_day") or (
                                rep_speed * n_reps
                            )
                            task_speeds[t_id] = (
                                f"{rep_speed:>4.1f} ns/d [Agg: {agg_speed:>5.1f} ({n_reps}x)]"
                            )
                        else:
                            task_speeds[t_id] = payload.get("speed", "--")

                        task_etas[t_id] = payload.get("eta", "--:--")
                        task_nodes[t_id] = payload.get("node") or "node"

                        clean_name = t_id.strip()
                        m_chunk_wrap = re.match(
                            r"^RUN_MD_CHUNK_\d+_(.+?)\s+(\(#\d+\))\s*\([^)]*\)$",
                            clean_name,
                        )
                        if m_chunk_wrap:
                            clean_name = (
                                f"MD: {m_chunk_wrap.group(1)} {m_chunk_wrap.group(2)}"
                            )
                        elif ":" in clean_name and not clean_name.startswith(
                            ("EQ:", "MD:", "SWISH:")
                        ):
                            clean_name = clean_name.split(":")[-1].strip()

                        if tool_raw not in ("mdrun", "gmx", "gromacs"):
                            tag_prefix = f"[{tool_raw.upper()}] "
                            if not clean_name.upper().startswith(f"[{tool_raw.upper()}]"):
                                clean_name = tag_prefix + clean_name

                        task_names[t_id] = clean_name

                    except (socket.timeout, BlockingIOError):
                        break
                    except Exception:
                        break

                for t_id in list(task_bars.keys()):
                    if now - last_seen.get(t_id, now) > 20.0:
                        progress.remove_task(task_bars[t_id])
                        del task_bars[t_id]
                        p_id = task_parents.pop(t_id, None)
                        if p_id and p_id in task_children:
                            task_children[p_id].discard(t_id)
                        for child_id in list(task_children.pop(t_id, set())):
                            if child_id in task_bars:
                                progress.remove_task(task_bars[child_id])
                                del task_bars[child_id]
                            task_parents.pop(child_id, None)
                            last_seen.pop(child_id, None)
                            task_real_step.pop(child_id, None)
                        last_seen.pop(t_id, None)
                        task_real_step.pop(t_id, None)

                for t_id, real_st in list(task_real_step.items()):
                    if now - last_seen.get(t_id, now) > 20.0:
                        continue

                    rate_s = task_rate_steps_s.get(t_id, 0.0)
                    time_delta = now - task_packet_time.get(t_id, now)
                    total_st = task_totals.get(t_id, 1)
                    dt_ps = task_dts.get(t_id)
                    is_sub = task_is_subtask.get(t_id, False)

                    if dt_ps is not None and dt_ps > 0:
                        if rate_s > 0 and time_delta > 0:
                            interpolated_st = min(
                                total_st, int(real_st + rate_s * time_delta)
                            )
                        else:
                            interpolated_st = real_st
                        t_str = f" [dim]({format_time_ns(interpolated_st * dt_ps)})[/dim]"
                        step_info = f"[cyan]{format_steps(interpolated_st)}/{format_steps(total_st)}[/cyan]{t_str}"
                    else:
                        interpolated_st = real_st
                        step_info = (
                            task_raw_step_info.get(t_id)
                            or f"[cyan]{format_steps(interpolated_st)}/{format_steps(total_st)}[/cyan]"
                        )

                    if is_sub:
                        lbl = task_subtask_label.get(t_id) or task_names.get(t_id, t_id)
                        desc = f"  [dim]└──[/dim] [bold cyan]{escape(lbl)}[/bold cyan]"
                        f_node_display = ""
                        f_gpu = ""
                        f_power = ""
                        f_speed = task_speeds.get(t_id, "")
                        f_eta = ""
                    else:
                        desc = f"[green]{escape(task_names.get(t_id, t_id))}[/green]"
                        f_node_display = f"({escape(task_nodes.get(t_id, 'node'))})"
                        f_gpu = escape(task_gpus.get(t_id, "GPU"))
                        f_power = task_powers.get(t_id, "-- W")
                        f_speed = task_speeds.get(t_id, "--")
                        f_eta = task_etas.get(t_id, "--:--")

                    if t_id not in task_bars:
                        task_bars[t_id] = progress.add_task(
                            desc,
                            total=total_st,
                            completed=interpolated_st,
                            node_display=f_node_display,
                            gpu=f_gpu,
                            power=f_power,
                            speed=f_speed,
                            step_info=step_info,
                            eta=f_eta,
                        )
                    else:
                        task_obj = progress._tasks[task_bars[t_id]]
                        if total_st is not None and task_obj.total != total_st:
                            progress.update(task_bars[t_id], total=total_st)

                        progress.update(
                            task_bars[t_id],
                            description=desc,
                            completed=interpolated_st,
                            node_display=f_node_display,
                            gpu=f_gpu,
                            power=f_power,
                            speed=f_speed,
                            step_info=step_info,
                            eta=f_eta,
                        )

                # Reorder progress._tasks so child tasks are strictly nested directly below their parent
                ordered_tasks = {}
                root_ids = [t for t in task_bars.keys() if not task_is_subtask.get(t, False)]
                for r_id in root_ids:
                    if r_id in task_bars:
                        r_bar_id = task_bars[r_id]
                        if r_bar_id in progress._tasks:
                            ordered_tasks[r_bar_id] = progress._tasks[r_bar_id]
                        for c_id in sorted(task_children.get(r_id, set())):
                            if c_id in task_bars:
                                c_bar_id = task_bars[c_id]
                                if c_bar_id in progress._tasks:
                                    ordered_tasks[c_bar_id] = progress._tasks[c_bar_id]
                for t_id, b_id in task_bars.items():
                    if b_id in progress._tasks and b_id not in ordered_tasks:
                        ordered_tasks[b_id] = progress._tasks[b_id]
                progress._tasks = ordered_tasks

                time.sleep(0.05)

    except KeyboardInterrupt:
        console.print("\n[yellow]Monitoring stopped.[/yellow]")
    finally:
        sock.close()
