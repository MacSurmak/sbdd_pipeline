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
    poses_name: str | None = None, top_n: int = 5, system_prefix: str = "hit"
) -> None:
    """
    =========================================================================================
    ARCHITECTURAL STUB NOTICE: 'sbdd select-hits' PENDING IMPLEMENTATION
    =========================================================================================
    CURRENT STATE:
      This command is currently a stub reporting pending status.
    REQUIRED CONTRACT:
      1. Query 'reports/ranking.csv' (or 'Pose.scores_dict' across candidate poses).
      2. Filter top N unique chemical ligands based on Ultimate_Score / GlideScore / MM-GBSA.
      3. For each selected pose:
         - Resolve its parent ReceptorFrame and SubstanceState.
         - Query existing cofactors on the parent frame.
         - Instantiate a new 3NF 'System' record in state 'raw' (ready for 'sbdd md-setup').
      4. Commit transaction to SQLite atomically without modifying pose or frame assets.
    =========================================================================================
    """
    session = init_db()
    console.print(
        f"[yellow]Notice: Full select-hits logic is pending implementation. Signature mapped successfully (poses={poses_name}, top={top_n}).[/yellow]"
    )
    session.close()
