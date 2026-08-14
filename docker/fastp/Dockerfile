FROM condaforge/miniforge3:latest

RUN conda install -y -c bioconda -c conda-forge fastp=0.23.* && \
    conda clean -afy

WORKDIR /data

ENTRYPOINT ["fastp"]