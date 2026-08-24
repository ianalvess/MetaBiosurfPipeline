FROM condaforge/miniforge3:latest

RUN conda install -y -n base -c conda-forge mamba && \
    conda clean -afy

RUN mamba create -y -n gtdbtk -c bioconda -c conda-forge "gtdbtk>=2.7" && \
    mamba clean -afy

ENV PATH=/opt/conda/envs/gtdbtk/bin:$PATH

WORKDIR /data

ENTRYPOINT ["gtdbtk"]