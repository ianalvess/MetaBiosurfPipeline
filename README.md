# MetaBiosurfPipeline 

Automated metagenomics pipeline for taxonomic and functional identification of biosurfactant-producing organisms in extreme environments.

## Overview

This pipeline processes shotgun metagenomic sequencing data (Illumina short-read) from extreme environments (hypersaline, hydrothermal, acidic/alkaline, etc.) through quality control, assembly, binning, quality assessment, taxonomic classification, and biosurfactant-focused functional annotation.

The pipeline is fully containerized: each tool runs in its own Docker image, orchestrated by a single Python CLI (`pipeline/cli.py`).

## Pipeline Steps

1. **QC** — [fastp](https://github.com/OpenGene/fastp): read trimming and quality filtering. Supports single-end and paired-end reads.
2. **Assembly** — [MEGAHIT](https://github.com/voutcn/megahit): de novo metagenome assembly.
3. **Mapping** — [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) + [samtools](https://github.com/samtools/samtools): maps clean reads back to assembled contigs for coverage estimation.
4. **Binning** — three binners run independently and are reconciled:
   - [MetaBAT2](https://bitbucket.org/berkeleylab/metabat)
   - [MaxBin2](https://sourceforge.net/projects/maxbin2/)
   - [CONCOCT](https://github.com/BinPro/CONCOCT)
   - [DAS Tool](https://github.com/cmks/DAS_Tool): reconciles the three binner outputs into a single non-redundant set of high-confidence bins.
5. **CheckM2** — evaluates completeness and contamination of each final bin using a machine-learning model (more robust than marker-gene counting, especially for divergent extremophile lineages).
6. **Kraken2** — whole-sample taxonomic classification, run directly on the clean reads (not per-bin). Uses a memory-capped database (`k2_standard_16gb`).
7. **BioSurfDB search** — per-bin gene prediction (Prodigal) followed by a DIAMOND search against a local [BioSurfDB](https://www.biosurfdb.org/) database, mapping hits to biosurfactant biosynthesis pathway categories (surfactin, rhamnolipids, lipopeptides, etc.).

## Requirements

- Docker Desktop (with WSL2 backend on Windows)
- Python 3.11 (conda environment `metagen-biosurf`)
- ~30GB+ free disk space for reference databases
- 16GB+ RAM recommended (32GB+ preferred; Kraken2's capped database was specifically chosen to work within this constraint)

### Python environment

```bash
conda env create -f environment.yml
conda activate metagen-biosurf
```

## Project Structure

```
metagen-biosurf/
├── config/
│   └── example.yaml          # per-sample configuration
├── data/
│   ├── <raw fastq files>
│   ├── checkm2_db/           # CheckM2 reference database
│   ├── kraken2_db/           # Kraken2 capped database (k2_standard_16gb)
│   └── biosurfdb/            # BioSurfDB DIAMOND database + mapping files
├── docker/
│   ├── fastp/
│   ├── megahit/
│   ├── mapping/              # bwa-mem2 + samtools
│   ├── metabat2/
│   ├── maxbin2/
│   ├── concoct/
│   ├── dastool/              # also provides prodigal + diamond, reused by BioSurfDB step
│   ├── checkm2/
│   └── kraken2/
├── pipeline/
│   └── cli.py                # orchestrator
├── results/
│   └── <sample_id>/          # all outputs, one folder per sample
└── environment.yml
```

## Configuration

Each sample is described by a YAML file (see `config/example.yaml`):

```yaml
sample_id: SRR328983
environment_type: hypersaline

reads:
  r1: data/SRR328983.fastq.gz
  # r2: data/SRR328983_R2.fastq.gz   # omit for single-end

output_dir: results/SRR328983

steps:
  qc:
    docker_image: metagen-biosurf/fastp
  assembly:
    docker_image: metagen-biosurf/megahit
  mapping:
    docker_image: metagen-biosurf/mapping
  binning_metabat2:
    docker_image: metagen-biosurf/metabat2
  binning_maxbin2:
    docker_image: metagen-biosurf/maxbin2
  binning_concoct:
    docker_image: metagen-biosurf/concoct
  binning_reconcile:
    docker_image: metagen-biosurf/dastool
  checkm2:
    docker_image: metagen-biosurf/checkm2
    db_path: data/checkm2_db
  kraken2:
    docker_image: metagen-biosurf/kraken2
    db_path: data/kraken2_db
  biosurfdb:
    docker_image: metagen-biosurf/dastool
    db_path: data/biosurfdb
```

## Usage

### 1. Build the Docker images (one-time, or after a Dockerfile change)

```bash
docker build -t metagen-biosurf/fastp docker/fastp
docker build -t metagen-biosurf/megahit docker/megahit
docker build -t metagen-biosurf/mapping docker/mapping
docker build -t metagen-biosurf/metabat2 docker/metabat2
docker build -t metagen-biosurf/maxbin2 docker/maxbin2
docker build -t metagen-biosurf/concoct docker/concoct
docker build -t metagen-biosurf/dastool docker/dastool
docker build -t metagen-biosurf/checkm2 docker/checkm2
docker build -t metagen-biosurf/kraken2 docker/kraken2
```

### 2. Download reference databases (one-time)

**CheckM2:**
```bash
docker run --rm -v ${PWD}/data/checkm2_db:/db metagen-biosurf/checkm2 database --download --path /db
```

**Kraken2** (`k2_standard_16gb`):
```bash
mkdir data/kraken2_db
docker run --rm -v ${PWD}/data/kraken2_db:/db --entrypoint bash metagen-biosurf/kraken2 \
  -c "cd /db && wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16_GB_<date>.tar.gz && tar xvzf k2_standard_16_GB_<date>.tar.gz && rm k2_standard_16_GB_<date>.tar.gz"
```
(check the current filename at https://benlangmead.github.io/aws-indexes/k2)

**BioSurfDB:** no public bulk download exists. Provide your own local copy (DIAMOND-formatted `.dmnd` file plus `acc2biosurfdb.map` and `biosurfdb.map` mapping files) in `data/biosurfdb/`.

### 3. Run the pipeline

```bash
python pipeline/cli.py --config config/example.yaml
```

## Outputs

All outputs are written under `results/<sample_id>/`:

- `<sample>.clean.fastq.gz`, `fastp.json`, `fastp.html` — QC results
- `assembly/final.contigs.fa` — assembled contigs
- `mapping/sorted.bam` — read alignment to contigs
- `binning/dastool/dastool_DASTool_bins/*.fa` — final reconciled bins
- `checkm2/quality_report.tsv` — completeness/contamination per bin
- `kraken2/kraken2_report.txt` — full Kraken2 report
- `kraken2/report/top10_species.{csv,png}`, `top10_genus.{csv,png}` — top-10 taxa charts
- `functional/biosurfdb/summary.tsv` — all DIAMOND hits, mapped to pathway categories
- `functional/biosurfdb/report/top20_categories.{csv,png}` — top-20 categories, split by binner/identity confidence
- `functional/biosurfdb/report/global_sample_categories.{csv,png}` — top-20 categories, whole sample combined

## Notes and Known Issues

- **DAS Tool (docopt R package bug):** requires pinning `r-base=4.2` and `r-docopt=0.7.1` in `docker/dastool/Dockerfile` to avoid a known R reference-class incompatibility.
- **DAS Tool + CONCOCT contig ID mismatches:** CONCOCT's contig cut/merge process can occasionally produce contig IDs that don't exactly match the assembly FASTA. The pipeline sanitizes `scaffolds2bin.tsv` files against the assembly before calling DAS Tool.
- **MEGAHIT `--force`:** used to avoid Windows/WSL2 file-lock issues when an output directory from a previous run still exists.
- **GTDB-Tk was removed** from this pipeline after repeated out-of-memory failures (the `bac120` reference tree requires ~50GB RAM at the `pplacer` step, exceeding available 32GB). Kraken2 with a capped 16GB database was adopted instead, trading fine-grained per-bin taxonomy for a reliable whole-sample overview.
- **Docker/WSL2 disk space:** the Docker Desktop virtual disk (`ext4.vhdx`) can grow very large with big reference databases. If disk space runs low, consider moving the Docker Desktop disk image location to a drive with more free space (Settings → Resources → Advanced → Disk image location).
- **Transient Docker/WSL2 memory issues:** if a container fails with a broken pipe error (`GetOverlappedResult`) or a step that previously worked suddenly runs out of memory, try `wsl --shutdown` followed by restarting Docker Desktop before re-running — this has resolved several otherwise-unexplained failures during development.