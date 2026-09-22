"""
Parses and validates the YAML configuration file.
"""
import yaml
from pathlib import Path
from loguru import logger

def load_config(config_path: Path) -> dict:
    """
    Loads YAML configuration and returns a dictionary.
    """
    if not config_path.exists():
        logger.critical(f"Configuration file not found: {config_path}")
        raise FileNotFoundError(f"Missing {config_path}")

    with open(config_path, "r") as f:
        try:
            config = yaml.safe_load(f)
            logger.debug(f"Successfully loaded configuration from {config_path.name}")
            return config
        except yaml.YAMLError as e:
            logger.critical(f"Invalid YAML format in {config_path.name}: {e}")
            raise