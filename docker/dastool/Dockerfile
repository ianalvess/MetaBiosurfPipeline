FROM condaforge/miniforge3:latest

RUN conda create -y -n dastool -c bioconda -c conda-forge das_tool && \
    conda clean -afy

ENV PATH=/opt/conda/envs/dastool/bin:$PATH

WORKDIR /data

ENTRYPOINT ["DAS_Tool"]