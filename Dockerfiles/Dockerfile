FROM continuumio/miniconda3
ENV TOOL=fastcov
ENV VERSION=0.1.3

RUN apt update && apt install -y procps wget gzip && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda update -n base -c defaults conda
RUN conda config --add channels default && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge
 

RUN conda install mamba
RUN mamba install matplotlib pandas openpyxl pysam seaborn
RUN conda clean -a

RUN git clone -b $VERSION https://github.com/RaverJay/$TOOL.git
ENV PATH /fastcov/:$PATH

RUN fastcov.py -h
