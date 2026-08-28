FROM condaforge/miniforge3:latest

RUN conda install -y -n base -c conda-forge mamba && \
    conda clean -afy

RUN mamba create -y -n kraken2 -c bioconda -c conda-forge kraken2 krona && \
    mamba clean -afy

ENV PATH=/opt/conda/envs/kraken2/bin:$PATH

WORKDIR /data

ENTRYPOINT ["kraken2"]