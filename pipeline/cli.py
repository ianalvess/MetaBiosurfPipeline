"""Orchestrator for the metagenomics pipeline (biosurfactants / extreme environments).

Currently implements:
  1. config loading;
  2. Docker daemon check;
  3. QC step (fastp), run as a container. Supports single-end and paired-end.
  4. assembly step (MEGAHIT), run as a container. Supports single-end and paired-end.

Remaining steps (binning, taxonomy, functional annotation)
will be added one by one, each invoking its own container.
"""

import shutil
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


def is_paired(config: dict) -> bool:
    return "r2" in config["reads"] and config["reads"]["r2"]


def run_qc(client: docker.DockerClient, config: dict, project_root: Path) -> dict:
    qc_cfg = config["steps"]["qc"]
    image = qc_cfg["docker_image"]

    output_dir = project_root / config["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)

    paired = is_paired(config)
    r1 = config["reads"]["r1"]
    clean_r1_name = f"{Path(r1).stem}.clean.fastq.gz"

    container_cmd = [
        "-i", f"/data/{Path(r1).name}",
        "-o", f"/output/{clean_r1_name}",
        "-j", "/output/fastp.json",
        "-h", "/output/fastp.html",
    ]

    clean_reads = {"r1": output_dir / clean_r1_name}

    if paired:
        r2 = config["reads"]["r2"]
        clean_r2_name = f"{Path(r2).stem}.clean.fastq.gz"
        container_cmd += [
            "-I", f"/data/{Path(r2).name}",
            "-O", f"/output/{clean_r2_name}",
        ]
        clean_reads["r2"] = output_dir / clean_r2_name

    logger.info(
        f"Running QC step (fastp) with image '{image}' "
        f"({'paired-end' if paired else 'single-end'})..."
    )

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

    return clean_reads


def run_assembly(
    client: docker.DockerClient,
    config: dict,
    project_root: Path,
    clean_reads: dict,
) -> None:
    assembly_cfg = config["steps"]["assembly"]
    image = assembly_cfg["docker_image"]

    output_dir = project_root / config["output_dir"]
    assembly_output = output_dir / "assembly"

    # MEGAHIT refuses to run if the output directory already exists.
    if assembly_output.exists():
        logger.info(f"Removing existing assembly output at {assembly_output}...")
        shutil.rmtree(assembly_output)

    paired = "r2" in clean_reads

    if paired:
        container_cmd = [
            "-1", f"/output/{clean_reads['r1'].name}",
            "-2", f"/output/{clean_reads['r2'].name}",
            "-o", "/output/assembly",
        ]
    else:
        container_cmd = [
            "-r", f"/output/{clean_reads['r1'].name}",
            "-o", "/output/assembly",
        ]

    logger.info(
        f"Running assembly step (MEGAHIT) with image '{image}' "
        f"({'paired-end' if paired else 'single-end'})..."
    )

    volumes = {
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
    logger.success(f"Assembly step finished. Output in {assembly_output}")


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

    clean_reads = run_qc(client, config, project_root)
    run_assembly(client, config, project_root, clean_reads)

    # TODO: chain the remaining steps defined in config["steps"]


if __name__ == "__main__":
    main()