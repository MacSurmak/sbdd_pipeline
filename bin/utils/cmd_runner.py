"""
File: src/utils/cmd_runner.py
Description: Generic, tool-agnostic command execution engine with live log monitoring.
Uses modular parsers from src/core/parsers.py and outputs real-time parser_debug.log.
"""
import subprocess
import time
import shutil
from pathlib import Path
from loguru import logger

from core.parsers import get_parser, BaseLogParser

def run_command(cmd: list, cwd: Path, log_name: str, progress=None, task_id=None, 
                parser: BaseLogParser = None, stdin: Path = None, total_override: int = None,
                step_offset: int = 0) -> bool:
    cmd_str = " ".join(str(x) for x in cmd)
    
    log_dir = cwd / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{log_name}.log"
    
    logger.info(f"Executing: {cmd[0]} (Log: {log_file.name})")
    logger.debug(f"Full command: {cmd_str}")
    
    # Auto-resolve parser if not explicitly passed
    active_parser = parser or get_parser(cmd)
    
    # Target log detection: external file for Glide (.in), Prime (-jobname), or GROMACS (-deffnm)
    if len(cmd) > 1 and str(cmd[1]).endswith(".in"):
        target_job_log = cwd / f"{Path(cmd[1]).stem}.log"
        target_job_log.unlink(missing_ok=True)
    elif "-jobname" in cmd:
        job_name_arg = cmd[cmd.index("-jobname") + 1]
        target_job_log = cwd / f"{job_name_arg}.log"
    elif "-deffnm" in cmd:
        deffnm_arg = cmd[cmd.index("-deffnm") + 1]
        target_job_log = cwd / f"{deffnm_arg}.log"
    else:
        target_job_log = log_file

    start_time = time.time()

    try:
        with open(log_file, "w") as f_out:
            f_out.write(f"COMMAND: {cmd_str}\n\nOUTPUT:\n")
            f_out.flush()

            f_in = open(stdin, "r") if stdin and stdin.exists() else None
            process = subprocess.Popen(
                cmd, cwd=cwd, stdin=f_in, stdout=f_out, stderr=subprocess.STDOUT
            )

            while process.poll() is None:
                # Dynamically resolve active log file (handles GROMACS .part000X.log suffixes)
                active_log = target_job_log
                if "-deffnm" in cmd:
                    deffnm_arg = cmd[cmd.index("-deffnm") + 1]
                    matching_logs = sorted(cwd.glob(f"{deffnm_arg}*.log"), key=lambda p: p.stat().st_mtime, reverse=True)
                    if matching_logs:
                        active_log = matching_logs[0]

                if progress and (task_id is not None) and active_parser and active_log.exists():
                    try:
                        with open(active_log, "rb") as f_read:
                            f_read.seek(0, 2)
                            f_size = f_read.tell()
                            f_read.seek(max(0, f_size - 8192))
                            raw_data = f_read.read().decode("utf-8", errors="ignore")

                        update = active_parser.parse(raw_data)
                        if update:
                            update_kwargs = {}
                            if update.completed is not None:
                                current_step = max(0, update.completed - step_offset)
                                update_kwargs["completed"] = current_step
                            
                            if total_override is not None:
                                update_kwargs["total"] = total_override
                            elif update.total is not None: 
                                update_kwargs["total"] = update.total
                                
                            if update.description is not None: update_kwargs["description"] = update.description
                            
                            progress.update(task_id, **update_kwargs)
                    except Exception:
                        pass

                time.sleep(0.5)

        returncode = process.returncode
        
        # Explicit progress finalization
        if progress and (task_id is not None):
            task = progress._tasks.get(task_id)
            if returncode == 0 and task and task.total:
                progress.update(task_id, completed=task.total)
            elif returncode != 0:
                progress.update(task_id, description=f"[red]FAILED: {Path(cmd[0]).name}[/red]")

        # Move stray logs created by Schrodinger JobControl
        for stray_log in cwd.glob("*.log"):
            if stray_log.is_file() and stray_log.name != log_file.name:
                shutil.move(str(stray_log), str(log_dir / stray_log.name))

    except Exception as e:
        logger.critical(f"Failed to launch command {cmd[0]}: {e}")
        return False

    duration = time.time() - start_time
    
    if returncode != 0:
        logger.error(f"Command '{cmd[0]}' failed with exit code {returncode}.")
        try:
            with open(log_file, "r") as f:
                logger.error(f"Tail of output:\n{''.join(f.readlines()[-10:])}")
        except Exception:
            pass
        return False
        
    logger.success(f"Command '{cmd[0]}' completed in {duration:.1f}s.")
    return True