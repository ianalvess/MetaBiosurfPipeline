FROM condaforge/miniforge3:latest

RUN conda create -y -n mapping -c bioconda -c conda-forge bwa-mem2 samtools && \
    conda clean -afy

ENV PATH=/opt/conda/envs/mapping/bin:$PATH

WORKDIR /data

ENTRYPOINT ["/bin/bash"]