FROM condaforge/miniforge3:latest

RUN conda create -y -n concoct -c bioconda -c conda-forge concoct && \
    conda clean -afy

ENV PATH=/opt/conda/envs/concoct/bin:$PATH

WORKDIR /data

ENTRYPOINT ["/bin/bash"]