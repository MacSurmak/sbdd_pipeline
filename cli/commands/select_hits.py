"""
File: cli/commands/select_hits.py
Description: CLI handler to select top candidate hits from docking and promote them to Systems for MD.
"""

import sys
from pathlib import Path
from rich.console import Console
from cli.core.db import init_db

console = Console()


def run_select_hits(
    poses_name: str | None = None,
    top_n: int = 5,
    system_prefix: str = "hit"
) -> None:
    """Selects top candidate hits from docking rankings and promotes to Systems."""
    session = init_db()
    console.print(f"[yellow]Notice: Full select-hits logic is pending implementation. Signature mapped successfully (poses={poses_name}, top={top_n}).[/yellow]")
    session.close()
