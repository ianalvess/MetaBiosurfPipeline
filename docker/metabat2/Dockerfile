FROM condaforge/miniforge3:latest

RUN conda create -y -n metabat2 -c bioconda -c conda-forge metabat2 && \
    conda clean -afy

ENV PATH=/opt/conda/envs/metabat2/bin:$PATH

WORKDIR /data

ENTRYPOINT ["metabat2"]