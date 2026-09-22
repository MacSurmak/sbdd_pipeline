"""
Configures the Loguru logger for both console output (colorized) and file output.
"""
import sys
from pathlib import Path
from loguru import logger
from rich.console import Console

console = Console()

def setup_logger(log_file: Path = None, verbose: bool = False):
    """
    Initializes logger settings.
    """
    logger.remove()  # Remove default handler
    
    # Console handler
    log_level = "DEBUG" if verbose else "INFO"
    logger.add(
        lambda msg: console.print(msg, end=""),
        level=log_level,
        format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>"
    )
    
    # File handler (if working directory is initialized)
    if log_file:
        logger.add(
            str(log_file),
            level="DEBUG",
            format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {name}:{function}:{line} - {message}",
            rotation="10 MB"
        )
