"""Orchestrator for the metagenomics pipeline (biosurfactants / extreme environments).

For now, this module only validates that:
  1. the configuration file is read correctly;
  2. the Docker daemon is reachable.

The actual steps (QC, assembly, binning, taxonomy, functional annotation)
will be added one by one, each invoking its own container.
"""

import sys
from pathlib import Path

import click
import docker
import yaml
from loguru import logger


def load_config(config_path: Path) -> dict:
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def check_docker() -> bool:
    try:
        client = docker.from_env()
        client.ping()
        return True
    except Exception as exc:  # noqa: BLE001
        logger.error(f"Could not reach the Docker daemon: {exc}")
        return False


@click.command()
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to the sample's YAML configuration file.",
)
def main(config_path: Path) -> None:
    logger.info(f"Loading configuration from: {config_path}")
    config = load_config(config_path)
    logger.info(f"Sample: {config.get('sample_id')} | Environment: {config.get('environment_type')}")

    logger.info("Checking Docker access...")
    if not check_docker():
        sys.exit(1)
    logger.success("Docker is reachable. Environment ready to run the pipeline.")

    # TODO: chain the steps defined in config["steps"]


if __name__ == "__main__":
    main()