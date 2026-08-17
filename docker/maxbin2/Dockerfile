FROM condaforge/miniforge3:latest

RUN conda create -y -n maxbin2 -c bioconda -c conda-forge maxbin2 && \
    conda clean -afy

ENV PATH=/opt/conda/envs/maxbin2/bin:$PATH

WORKDIR /data

ENTRYPOINT ["run_MaxBin.pl"]