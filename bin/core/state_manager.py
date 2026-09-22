"""
Validates dependencies and checkpoint files to ensure modules run in correct order.
"""
import sys
from pathlib import Path
from loguru import logger

class StateManager:
    @staticmethod
    def require_file(filepath: Path, module_name: str) -> None:
        """
        Validates existence of a required file. Exits gracefully if missing.
        """
        if not filepath.exists():
            logger.critical(f"Missing required dependency: {filepath}")
            logger.error(f"Cannot execute '{module_name}'. Please run prerequisites first.")
            sys.exit(1)

    @staticmethod
    def check_output_exists(filepath: Path, force: bool = False) -> bool:
        """
        Checks if output already exists to handle checkpoints / overwrites.
        """
        if filepath.exists():
            if force:
                logger.warning(f"Forcing overwrite of existing output: {filepath}")
                return False
            else:
                logger.info(f"Checkpoint found: {filepath.name} exists. Skipping step. Use --force to override.")
                return True
        return False