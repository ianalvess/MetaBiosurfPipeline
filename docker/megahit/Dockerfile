FROM condaforge/miniforge3:latest

RUN conda install -y -c bioconda -c conda-forge megahit=1.2.* && \
    conda clean -afy

WORKDIR /data

ENTRYPOINT ["megahit"]