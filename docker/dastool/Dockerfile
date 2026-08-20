FROM condaforge/miniforge3:latest

RUN conda install -y -n base -c conda-forge mamba && \
    conda clean -afy

RUN mamba create -y -n dastool -c bioconda -c conda-forge \
    "das_tool=1.1.6" \
    "r-base=4.2" \
    "r-docopt=0.7.1" && \
    mamba clean -afy

ENV PATH=/opt/conda/envs/dastool/bin:$PATH

WORKDIR /data

ENTRYPOINT ["DAS_Tool"]