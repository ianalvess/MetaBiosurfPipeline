FROM condaforge/miniforge3:latest

RUN conda install -y -n base -c conda-forge mamba && \
    conda clean -afy

RUN mamba create -y -n antismash -c bioconda -c conda-forge antismash && \
    mamba clean -afy

ENV PATH=/opt/conda/envs/antismash/bin:$PATH

WORKDIR /data

ENTRYPOINT ["antismash"]