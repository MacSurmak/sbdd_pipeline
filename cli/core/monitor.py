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
                    dt_ps = msg.get("dt_ps")
                    start_step = msg.get("start")
                    total_steps = msg.get("total")
                    current_step = msg.get("step")

                    if current_step is not None:
                        # Distinguish discrete pipeline chunks from continuous/external runs
                        is_chunk = bool(re.search(r"(?:_C\d+|RUN_MD_CHUNK)", task_name))

                        if (
                            is_chunk
                            and start_step is not None
                            and total_steps is not None
                            and total_steps > start_step
                        ):
                            display_total = total_steps - start_step
                            display_step = max(
                                0, min(display_total, current_step - start_step)
                            )
                        else:
                            display_total = total_steps
                            display_step = current_step

                        # Physical simulation time string (absolute trajectory time: step * dt)
                        time_info = ""
                        if dt_ps is not None and current_step is not None:
                            t_str = format_time_ns(current_step * dt_ps)
                            if t_str:
                                time_info = f" [dim]({t_str})[/dim]"

                        if display_total is not None and display_step is not None:
                            step_info_str = f"[cyan]{format_steps(display_step)}/{format_steps(display_total)}[/cyan]{time_info}"
                        elif display_step is not None:
                            step_info_str = (
                                f"[cyan]{format_steps(display_step)}[/cyan]{time_info}"
                            )
                        else:
                            step_info_str = time_info

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

                        rem_steps = (
                            (display_total - display_step)
                            if (display_total is not None and display_step is not None)
                            else None
                        )
                        eta_display = "--:--"
                        if rem_steps is not None and rem_steps > 0:
                            ema = task_stats[task_name]["ema_rate"]
                            if ema and ema > 0:
                                eta_display = format_eta(rem_steps / ema)
                            else:
                                eta_display = "calculating..."
                        elif rem_steps == 0:
                            eta_display = "done"

                        if task_name not in task_bars:
                            task_bars[task_name] = progress.add_task(
                                f"[green]{escape(clean_name)}",
                                total=display_total,
                                completed=display_step,
                                node=escape(node_name),
                                speed=speed_display,
                                step_info=step_info_str,
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
                                step_info=step_info_str,
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
