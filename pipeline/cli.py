"""Orchestrator for the metagenomics pipeline (biosurfactants / extreme environments).

Currently implements:
  1. config loading;
  2. Docker daemon check;
  3. QC step (fastp), run as a container. Supports single-end and paired-end.
  4. assembly step (MEGAHIT), run as a container. Supports single-end and paired-end.
  5. mapping step (bwa-mem2 + samtools): maps clean reads back to contigs.
  6. binning step: MetaBAT2, MaxBin2, and CONCOCT run independently, then
     reconciled with DAS Tool into a final, non-redundant set of bins.

Remaining steps (taxonomy, functional annotation)
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
) -> Path:
    assembly_cfg = config["steps"]["assembly"]
    image = assembly_cfg["docker_image"]

    output_dir = project_root / config["output_dir"]
    assembly_output = output_dir / "assembly"

    paired = "r2" in clean_reads

    # --force lets MEGAHIT overwrite an existing output directory,
    # avoiding Windows/WSL2 file-lock issues when deleting it from Python.
    if paired:
        container_cmd = [
            "-1", f"/output/{clean_reads['r1'].name}",
            "-2", f"/output/{clean_reads['r2'].name}",
            "-o", "/output/assembly",
            "--force",
        ]
    else:
        container_cmd = [
            "-r", f"/output/{clean_reads['r1'].name}",
            "-o", "/output/assembly",
            "--force",
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

    return assembly_output / "final.contigs.fa"


def run_mapping(
    client: docker.DockerClient,
    config: dict,
    project_root: Path,
    clean_reads: dict,
) -> Path:
    mapping_cfg = config["steps"]["mapping"]
    image = mapping_cfg["docker_image"]

    output_dir = project_root / config["output_dir"]
    mapping_output = output_dir / "mapping"
    mapping_output.mkdir(parents=True, exist_ok=True)

    contigs = "/output/assembly/final.contigs.fa"
    sorted_bam = "/output/mapping/sorted.bam"

    reads_arg = f"/output/{clean_reads['r1'].name}"
    if "r2" in clean_reads:
        reads_arg += f" /output/{clean_reads['r2'].name}"

    script = (
        f"set -e && "
        f"bwa-mem2 index {contigs} && "
        f"bwa-mem2 mem -t 4 {contigs} {reads_arg} | "
        f"samtools sort -o {sorted_bam} - && "
        f"samtools index {sorted_bam}"
    )

    logger.info(f"Running mapping step (bwa-mem2 + samtools) with image '{image}'...")

    volumes = {
        str(output_dir.resolve()): {"bind": "/output", "mode": "rw"},
    }

    logs = client.containers.run(
        image,
        command=["-c", script],
        volumes=volumes,
        remove=True,
        stdout=True,
        stderr=True,
    )
    logger.info(logs.decode("utf-8", errors="replace"))
    logger.success(f"Mapping step finished. Output in {mapping_output}")

    return mapping_output / "sorted.bam"


def _run_bash(
    client: docker.DockerClient,
    image: str,
    script: str,
    volumes: dict,
    step_label: str,
) -> None:
    logger.info(f"Running {step_label} with image '{image}'...")
    logs = client.containers.run(
        image,
        command=["-c", script],
        entrypoint="/bin/bash",
        volumes=volumes,
        remove=True,
        stdout=True,
        stderr=True,
    )
    logger.info(logs.decode("utf-8", errors="replace"))
    logger.success(f"{step_label} finished.")


def run_binning_metabat2(
    client: docker.DockerClient,
    config: dict,
    project_root: Path,
) -> None:
    binning_cfg = config["steps"]["binning_metabat2"]
    image = binning_cfg["docker_image"]

    output_dir = project_root / config["output_dir"]
    (output_dir / "binning" / "metabat2").mkdir(parents=True, exist_ok=True)

    contigs = "/output/assembly/final.contigs.fa"
    sorted_bam = "/output/mapping/sorted.bam"
    depth_file = "/output/mapping/depth.txt"

    volumes = {str(output_dir.resolve()): {"bind": "/output", "mode": "rw"}}

    logger.info("Computing contig depth for MetaBAT2/MaxBin2...")
    depth_logs = client.containers.run(
        image,
        command=["--outputDepth", depth_file, sorted_bam],
        entrypoint="jgi_summarize_bam_contig_depths",
        volumes=volumes,
        remove=True,
        stdout=True,
        stderr=True,
    )
    logger.info(depth_logs.decode("utf-8", errors="replace"))

    logger.info(f"Running MetaBAT2 binning with image '{image}'...")
    logs = client.containers.run(
        image,
        command=[
            "-i", contigs,
            "-a", depth_file,
            "-o", "/output/binning/metabat2/bin",
        ],
        volumes=volumes,
        remove=True,
        stdout=True,
        stderr=True,
    )
    logger.info(logs.decode("utf-8", errors="replace"))
    logger.success("MetaBAT2 binning finished.")


def run_binning_maxbin2(
    client: docker.DockerClient,
    config: dict,
    project_root: Path,
) -> None:
    binning_cfg = config["steps"]["binning_maxbin2"]
    image = binning_cfg["docker_image"]

    output_dir = project_root / config["output_dir"]
    (output_dir / "binning" / "maxbin2").mkdir(parents=True, exist_ok=True)

    volumes = {str(output_dir.resolve()): {"bind": "/output", "mode": "rw"}}

    # MaxBin2 wants a two-column abundance file: contigName  avgDepth.
    # Derive it from the depth.txt already produced for MetaBAT2 (col 1 and 3).
    script = (
        "set -e && "
        "tail -n +2 /output/mapping/depth.txt | "
        "awk '{print $1\"\\t\"$3}' > /output/binning/maxbin2/abundance.txt && "
        "run_MaxBin.pl "
        "-contig /output/assembly/final.contigs.fa "
        "-abund /output/binning/maxbin2/abundance.txt "
        "-out /output/binning/maxbin2/bin "
        "-thread 4"
    )

    _run_bash(client, image, script, volumes, "MaxBin2 binning")


def run_binning_concoct(
    client: docker.DockerClient,
    config: dict,
    project_root: Path,
) -> None:
    binning_cfg = config["steps"]["binning_concoct"]
    image = binning_cfg["docker_image"]

    output_dir = project_root / config["output_dir"]
    (output_dir / "binning" / "concoct" / "fasta_bins").mkdir(parents=True, exist_ok=True)

    volumes = {str(output_dir.resolve()): {"bind": "/output", "mode": "rw"}}

    script = (
        "set -e && "
        "cut_up_fasta.py /output/assembly/final.contigs.fa -c 10000 -o 0 "
        "--merge_last -b /output/binning/concoct/contigs_10K.bed "
        "> /output/binning/concoct/contigs_10K.fa && "
        "concoct_coverage_table.py /output/binning/concoct/contigs_10K.bed "
        "/output/mapping/sorted.bam > /output/binning/concoct/coverage_table.tsv && "
        "concoct --composition_file /output/binning/concoct/contigs_10K.fa "
        "--coverage_file /output/binning/concoct/coverage_table.tsv "
        "-b /output/binning/concoct/ && "
        "merge_cutup_clustering.py /output/binning/concoct/clustering_gt1000.csv "
        "> /output/binning/concoct/clustering_merged.csv && "
        "extract_fasta_bins.py /output/assembly/final.contigs.fa "
        "/output/binning/concoct/clustering_merged.csv "
        "--output_path /output/binning/concoct/fasta_bins"
    )

    _run_bash(client, image, script, volumes, "CONCOCT binning")


def run_binning_reconcile(
    client: docker.DockerClient,
    config: dict,
    project_root: Path,
) -> None:
    dastool_cfg = config["steps"]["binning_reconcile"]
    image = dastool_cfg["docker_image"]

    output_dir = project_root / config["output_dir"]
    (output_dir / "binning" / "dastool").mkdir(parents=True, exist_ok=True)

    volumes = {str(output_dir.resolve()): {"bind": "/output", "mode": "rw"}}

    # Convert each binner's fasta bins into DAS Tool's contigs2bin.tsv format.
    # Note: the conversion script was renamed from Fasta_to_Scaffolds2Bin.sh
    # to Fasta_to_Contig2Bin.sh in DAS Tool v1.1.4+.
    #
    # After conversion, sanitize each scaffolds2bin.tsv: some binners
    # (notably CONCOCT, after cutting contigs and merging clusters back)
    # can emit contig IDs that don't exactly match the assembly fasta due
    # to remapping edge cases. DAS Tool fails hard on any mismatch, so we
    # drop unmatched lines before handing the tables to it.
    script = (
        "set -e && "
        "Fasta_to_Contig2Bin.sh -i /output/binning/metabat2 -e fa "
        "> /output/binning/metabat2.scaffolds2bin.tsv && "
        "Fasta_to_Contig2Bin.sh -i /output/binning/maxbin2 -e fasta "
        "> /output/binning/maxbin2.scaffolds2bin.tsv && "
        "Fasta_to_Contig2Bin.sh -i /output/binning/concoct/fasta_bins -e fa "
        "> /output/binning/concoct.scaffolds2bin.tsv && "
        "grep '^>' /output/assembly/final.contigs.fa | sed 's/^>//; s/ .*//' "
        "> /output/binning/valid_contigs.txt && "
        "for f in metabat2 maxbin2 concoct; do "
        "  awk 'NR==FNR{valid[$1]=1; next} $1 in valid' "
        "  /output/binning/valid_contigs.txt "
        "  /output/binning/${f}.scaffolds2bin.tsv "
        "  > /output/binning/${f}.scaffolds2bin.clean.tsv && "
        "  mv /output/binning/${f}.scaffolds2bin.clean.tsv "
        "  /output/binning/${f}.scaffolds2bin.tsv; "
        "done && "
        "DAS_Tool "
        "-i /output/binning/metabat2.scaffolds2bin.tsv,"
        "/output/binning/maxbin2.scaffolds2bin.tsv,"
        "/output/binning/concoct.scaffolds2bin.tsv "
        "-l metabat2,maxbin2,concoct "
        "-c /output/assembly/final.contigs.fa "
        "-o /output/binning/dastool/dastool "
        "--write_bins "
        "-t 4"
    )

    _run_bash(client, image, script, volumes, "DAS Tool bin reconciliation")


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
    run_mapping(client, config, project_root, clean_reads)

    run_binning_metabat2(client, config, project_root)
    run_binning_maxbin2(client, config, project_root)
    run_binning_concoct(client, config, project_root)
    run_binning_reconcile(client, config, project_root)

    # TODO: chain the remaining steps (CheckM2, GTDB-Tk, functional annotation)


if __name__ == "__main__":
    main()