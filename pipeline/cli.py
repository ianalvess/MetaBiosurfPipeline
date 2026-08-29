"""Orchestrator for the metagenomics pipeline (biosurfactants / extreme environments).

Currently implements:
  1. config loading;
  2. Docker daemon check;
  3. QC step (fastp), run as a container. Supports single-end and paired-end.
  4. assembly step (MEGAHIT), run as a container. Supports single-end and paired-end.
  5. mapping step (bwa-mem2 + samtools): maps clean reads back to contigs.
  6. binning step: MetaBAT2, MaxBin2, and CONCOCT run independently, then
     reconciled with DAS Tool into a final, non-redundant set of bins.
  7. CheckM2: evaluates completeness/contamination of the final bins.
  8. Kraken2: whole-sample taxonomic classification run directly on the
     clean reads (low-memory, fits comfortably in 32GB RAM), with a
     publication-ready graphical summary of the top taxa.
  9. BioSurfDB search: Prodigal gene prediction per bin + DIAMOND search
     against a local BioSurfDB database, summarized by biosurfactant
     pathway category, with two reports: a top-20 table/chart split by
     identity confidence, and a whole-sample (all bins combined) view.
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


def run_checkm2(
    client: docker.DockerClient,
    config: dict,
    project_root: Path,
) -> None:
    checkm2_cfg = config["steps"]["checkm2"]
    image = checkm2_cfg["docker_image"]
    db_path = project_root / checkm2_cfg["db_path"]

    output_dir = project_root / config["output_dir"]
    checkm2_output = output_dir / "checkm2"

    volumes = {
        str(output_dir.resolve()): {"bind": "/output", "mode": "rw"},
        str(db_path.resolve()): {"bind": "/db", "mode": "ro"},
    }

    logger.info(f"Running CheckM2 with image '{image}'...")

    container_cmd = [
        "predict",
        "--input", "/output/binning/dastool/dastool_DASTool_bins",
        "--output-directory", "/output/checkm2",
        "--database_path", "/db/CheckM2_database/uniref100.KO.1.dmnd",
        "-x", "fa",
        "--force",
        "-t", "4",
    ]

    logs = client.containers.run(
        image,
        command=container_cmd,
        volumes=volumes,
        remove=True,
        stdout=True,
        stderr=True,
    )
    logger.info(logs.decode("utf-8", errors="replace"))
    logger.success(f"CheckM2 finished. Output in {checkm2_output}")


def run_kraken2(
    client: docker.DockerClient,
    config: dict,
    project_root: Path,
    clean_reads: dict,
) -> None:
    """Run Kraken2 directly on the clean reads for a whole-sample taxonomic
    classification. Low memory footprint (capped database), no dependency
    on binning/assembly for this step.
    """
    kraken2_cfg = config["steps"]["kraken2"]
    image = kraken2_cfg["docker_image"]
    db_path = project_root / kraken2_cfg["db_path"]

    output_dir = project_root / config["output_dir"]
    kraken2_output = output_dir / "kraken2"
    kraken2_output.mkdir(parents=True, exist_ok=True)

    volumes = {
        str(output_dir.resolve()): {"bind": "/output", "mode": "rw"},
        str(db_path.resolve()): {"bind": "/db", "mode": "ro"},
    }

    container_cmd = [
        "--db", "/db",
        "--threads", "4",
        "--report", "/output/kraken2/kraken2_report.txt",
        "--output", "/output/kraken2/kraken2_output.txt",
        "--use-names",
    ]

    if "r2" in clean_reads:
        container_cmd += [
            "--paired",
            f"/output/{clean_reads['r1'].name}",
            f"/output/{clean_reads['r2'].name}",
        ]
    else:
        container_cmd += [f"/output/{clean_reads['r1'].name}"]

    logger.info(f"Running Kraken2 with image '{image}'...")

    logs = client.containers.run(
        image,
        command=container_cmd,
        volumes=volumes,
        remove=True,
        stdout=True,
        stderr=True,
    )
    logger.info(logs.decode("utf-8", errors="replace"))
    logger.success(f"Kraken2 finished. Output in {kraken2_output}")


def generate_kraken2_report(config: dict, project_root: Path) -> None:
    """Build publication-ready horizontal bar charts of the top species-
    and genus-level taxa found in the whole sample by Kraken2, plus CSV
    tables for each. Pure Python — no container needed.
    """
    import pandas as pd
    import matplotlib.pyplot as plt

    output_dir = project_root / config["output_dir"]
    report_path = output_dir / "kraken2" / "kraken2_report.txt"
    report_dir = output_dir / "kraken2" / "report"
    report_dir.mkdir(parents=True, exist_ok=True)

    if not report_path.exists():
        logger.warning(f"Kraken2 report not found at {report_path} — skipping taxonomy report.")
        return

    df = pd.read_csv(
        report_path,
        sep="\t",
        header=None,
        names=["pct", "frag_total", "frag_direct", "rank", "taxid", "name"],
    )
    df["name"] = df["name"].str.strip()

    def build_chart(rank_code: str, rank_label: str, filename_stem: str) -> None:
        subset = df[df["rank"] == rank_code].sort_values("pct", ascending=False).head(10)
        if subset.empty:
            logger.warning(f"No {rank_label}-level taxa found in Kraken2 report — skipping chart.")
            return

        table_path = report_dir / f"{filename_stem}.csv"
        subset[["name", "pct", "frag_total"]].rename(
            columns={"pct": "percentage_of_reads", "frag_total": "read_count"}
        ).to_csv(table_path, index=False)
        logger.success(f"Top-10 {rank_label} table written to {table_path}")

        fig, ax = plt.subplots(figsize=(9, 5.5))
        ax.barh(
            subset["name"][::-1],
            subset["pct"][::-1],
            color="#4575b4",
            edgecolor="white",
            linewidth=0.5,
        )
        ax.set_xlabel("Percentage of classified reads (%)", fontsize=11)
        ax.set_ylabel("")
        ax.set_title(
            f"Kraken2 Taxonomic Classification — Whole Sample\n"
            f"Top 10 {rank_label} by read percentage",
            fontsize=13,
            fontweight="bold",
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()

        chart_path = report_dir / f"{filename_stem}.png"
        fig.savefig(chart_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        logger.success(f"{rank_label.capitalize()} chart written to {chart_path}")

    build_chart("S", "species", "top10_species")
    build_chart("G", "genus", "top10_genus")


def run_biosurfdb_search(
    client: docker.DockerClient,
    config: dict,
    project_root: Path,
) -> None:
    """Predict genes per bin (Prodigal) and search them against the
    local BioSurfDB DIAMOND database for biosurfactant-related hits.
    Reuses the dastool image, which already has prodigal and diamond.
    """
    biosurfdb_cfg = config["steps"]["biosurfdb"]
    image = biosurfdb_cfg["docker_image"]
    db_path = project_root / biosurfdb_cfg["db_path"]

    output_dir = project_root / config["output_dir"]
    functional_dir = output_dir / "functional" / "biosurfdb"
    (functional_dir / "proteins").mkdir(parents=True, exist_ok=True)
    (functional_dir / "hits").mkdir(parents=True, exist_ok=True)

    volumes = {
        str(output_dir.resolve()): {"bind": "/output", "mode": "rw"},
        str(db_path.resolve()): {"bind": "/db", "mode": "ro"},
    }

    script = (
        "set -e && "
        "for f in /output/binning/dastool/dastool_DASTool_bins/*.fa; do "
        '  name=$(basename "$f" .fa); '
        '  prodigal -i "$f" -a /output/functional/biosurfdb/proteins/${name}.faa '
        "  -p single -q; "
        "  diamond blastp "
        "    -q /output/functional/biosurfdb/proteins/${name}.faa "
        "    -d /db/biosurfdb.dmnd "
        "    -o /output/functional/biosurfdb/hits/${name}.tsv "
        "    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle "
        "    --evalue 1e-5 --max-target-seqs 1 --threads 4; "
        "done"
    )

    _run_bash(client, image, script, volumes, "BioSurfDB gene prediction + DIAMOND search")


def summarize_biosurfdb_hits(config: dict, project_root: Path) -> None:
    """Parse DIAMOND hits against BioSurfDB and map them to biosurfactant
    pathway categories, producing a summary table. Pure Python — runs on
    the host, no container needed.
    """
    biosurfdb_cfg = config["steps"]["biosurfdb"]
    db_path = project_root / biosurfdb_cfg["db_path"]

    output_dir = project_root / config["output_dir"]
    hits_dir = output_dir / "functional" / "biosurfdb" / "hits"
    summary_path = output_dir / "functional" / "biosurfdb" / "summary.tsv"

    # accession -> category ID
    acc2id: dict[str, str] = {}
    with open(db_path / "acc2biosurfdb.map", "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                acc2id[parts[0]] = parts[1]

    # category ID -> category name
    id2name: dict[str, str] = {}
    with open(db_path / "biosurfdb.map", "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                id2name[parts[0]] = parts[1]

    rows = []
    for hits_file in sorted(hits_dir.glob("*.tsv")):
        bin_name = hits_file.stem
        with open(hits_file, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    continue
                qseqid, sseqid, pident, length, evalue, bitscore = fields[:6]
                category_id = acc2id.get(sseqid, "")
                category_name = id2name.get(category_id, "unknown")
                rows.append(
                    {
                        "bin": bin_name,
                        "query_gene": qseqid,
                        "subject_accession": sseqid,
                        "pident": pident,
                        "evalue": evalue,
                        "bitscore": bitscore,
                        "category_id": category_id,
                        "category_name": category_name,
                    }
                )

    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(
            "bin\tquery_gene\tsubject_accession\tpident\tevalue\tbitscore\t"
            "category_id\tcategory_name\n"
        )
        for row in rows:
            f.write(
                f"{row['bin']}\t{row['query_gene']}\t{row['subject_accession']}\t"
                f"{row['pident']}\t{row['evalue']}\t{row['bitscore']}\t"
                f"{row['category_id']}\t{row['category_name']}\n"
            )

    logger.success(f"BioSurfDB summary written to {summary_path} ({len(rows)} hits)")


def generate_biosurfdb_report(config: dict, project_root: Path) -> None:
    """Build a publication-ready summary table and stacked bar chart from
    the BioSurfDB hits: top 20 functional categories by hit count, with
    each bar split by DIAMOND identity confidence band.
    """
    import pandas as pd
    import matplotlib.pyplot as plt

    output_dir = project_root / config["output_dir"]
    summary_path = output_dir / "functional" / "biosurfdb" / "summary.tsv"
    report_dir = output_dir / "functional" / "biosurfdb" / "report"
    report_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(summary_path, sep="\t")
    df["pident"] = df["pident"].astype(float)

    def confidence_band(pident: float) -> str:
        if pident < 30:
            return "<30% identity"
        elif pident < 50:
            return "30-50% identity"
        elif pident < 70:
            return "50-70% identity"
        else:
            return "\u226570% identity"

    df["confidence"] = df["pident"].apply(confidence_band)

    total_hits = len(df)
    top20 = (
        df.groupby("category_name")
        .size()
        .sort_values(ascending=False)
        .head(20)
        .index
    )

    # Percentage table
    table = (
        df[df["category_name"].isin(top20)]
        .groupby("category_name")
        .size()
        .reindex(top20)
        .reset_index(name="hit_count")
    )
    table["percentage_of_total_hits"] = (
        table["hit_count"] / total_hits * 100
    ).round(2)
    table_path = report_dir / "top20_categories.csv"
    table.to_csv(table_path, index=False)
    logger.success(f"Top-20 category table written to {table_path}")

    # Stacked bar chart: category x confidence band, as % of total hits
    pivot = (
        df[df["category_name"].isin(top20)]
        .groupby(["category_name", "confidence"])
        .size()
        .unstack(fill_value=0)
        .reindex(top20)
    )
    pivot_pct = pivot / total_hits * 100

    band_order = ["<30% identity", "30-50% identity", "50-70% identity", "\u226570% identity"]
    band_colors = ["#d73027", "#fc8d59", "#91bfdb", "#4575b4"]
    pivot_pct = pivot_pct[[b for b in band_order if b in pivot_pct.columns]]

    fig, ax = plt.subplots(figsize=(10, 7))
    bottom = pd.Series(0.0, index=pivot_pct.index)
    for band, color in zip(band_order, band_colors):
        if band not in pivot_pct.columns:
            continue
        ax.barh(
            pivot_pct.index,
            pivot_pct[band],
            left=bottom,
            label=band,
            color=color,
            edgecolor="white",
            linewidth=0.5,
        )
        bottom += pivot_pct[band]

    ax.set_xlabel("Percentage of total DIAMOND hits (%)", fontsize=11)
    ax.set_ylabel("")
    ax.set_title(
        "Top 20 BioSurfDB Functional Categories\n"
        "by Hit Percentage and Identity Confidence",
        fontsize=13,
        fontweight="bold",
    )
    ax.invert_yaxis()
    ax.legend(
        title="DIAMOND identity",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    chart_path = report_dir / "top20_categories_stacked.png"
    fig.savefig(chart_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.success(f"Stacked bar chart written to {chart_path}")


def generate_biosurfdb_global_report(config: dict, project_root: Path) -> None:
    """Build a whole-sample view of the BioSurfDB results: all bins
    combined into a single global picture (no per-bin or per-confidence
    breakdown), showing the top 20 functional categories as a percentage
    of the total sample's hits.
    """
    import pandas as pd
    import matplotlib.pyplot as plt

    output_dir = project_root / config["output_dir"]
    summary_path = output_dir / "functional" / "biosurfdb" / "summary.tsv"
    report_dir = output_dir / "functional" / "biosurfdb" / "report"
    report_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(summary_path, sep="\t")
    total_hits = len(df)

    counts = df.groupby("category_name").size().sort_values(ascending=False)
    top20 = counts.head(20)
    top20_pct = (top20 / total_hits * 100).round(2)

    table = top20.reset_index(name="hit_count")
    table["percentage_of_total_hits"] = top20_pct.values
    table_path = report_dir / "global_sample_categories.csv"
    table.to_csv(table_path, index=False)
    logger.success(f"Global sample category table written to {table_path}")

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.barh(
        top20_pct.index[::-1],
        top20_pct.values[::-1],
        color="#4575b4",
        edgecolor="white",
        linewidth=0.5,
    )
    ax.set_xlabel("Percentage of total sample hits (%)", fontsize=11)
    ax.set_ylabel("")
    ax.set_title(
        "BioSurfDB Functional Categories — Whole Sample\n"
        "Top 20 categories, all bins combined",
        fontsize=13,
        fontweight="bold",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    chart_path = report_dir / "global_sample_categories.png"
    fig.savefig(chart_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.success(f"Global sample chart written to {chart_path}")


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

    run_checkm2(client, config, project_root)

    run_kraken2(client, config, project_root, clean_reads)
    generate_kraken2_report(config, project_root)

    run_biosurfdb_search(client, config, project_root)
    summarize_biosurfdb_hits(config, project_root)
    generate_biosurfdb_report(config, project_root)
    generate_biosurfdb_global_report(config, project_root)


if __name__ == "__main__":
    main()