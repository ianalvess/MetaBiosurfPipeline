"""Orchestrator for the metagenomics pipeline (biosurfactants / extreme environments).

Currently implements:
  1. config loading;
  2. Docker daemon check;
  3. QC step (fastp), run as a container.

Remaining steps (assembly, binning, taxonomy, functional annotation)
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


def check_docker() -> docker.DockerClient | None:
    try:
        client = docker.from_env()
        client.ping()
        return client
    except Exception as exc:  # noqa: BLE001
        logger.error(f"Could not reach the Docker daemon: {exc}")
        return None


def run_qc(client: docker.DockerClient, config: dict, project_root: Path) -> None:
    qc_cfg = config["steps"]["qc"]
    image = qc_cfg["docker_image"]

    output_dir = project_root / config["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)

    r1 = config["reads"]["r1"]
    r2 = config["reads"]["r2"]

    logger.info(f"Running QC step (fastp) with image '{image}'...")

    container_cmd = [
        "-i", f"/data/{Path(r1).name}",
        "-I", f"/data/{Path(r2).name}",
        "-o", f"/output/{Path(r1).stem}.clean.fastq.gz",
        "-O", f"/output/{Path(r2).stem}.clean.fastq.gz",
        "-j", "/output/fastp.json",
        "-h", "/output/fastp.html",
    ]

    volumes = {
        str((project_root / "data").resolve()): {"bind": "/data", "mode": "ro"},
        str(output_dir.resolve()): {"bind": "/output", "mode": "rw"},
    }

    logs = client.containers.run(
        image,
        command=container_cmd,
        volumes=volumes,
        remove=True,
        stdout=True,
        stderr=True,
    )
    logger.info(logs.decode("utf-8", errors="replace"))
    logger.success(f"QC step finished. Output in {output_dir}")


@click.command()
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to the sample's YAML configuration file.",
)
def main(config_path: Path) -> None:
    project_root = Path(__file__).resolve().parent.parent

    logger.info(f"Loading configuration from: {config_path}")
    config = load_config(config_path)
    logger.info(f"Sample: {config.get('sample_id')} | Environment: {config.get('environment_type')}")

    logger.info("Checking Docker access...")
    client = check_docker()
    if client is None:
        sys.exit(1)
    logger.success("Docker is reachable.")

    run_qc(client, config, project_root)

    # TODO: chain the remaining steps defined in config["steps"]


if __name__ == "__main__":
    main()