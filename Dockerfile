# Set the base image to debian based miniconda2
FROM conda/miniconda3

# File Author / Maintainer
MAINTAINER JAY-Doq

SHELL ["/bin/bash","-c"]

COPY environments /environments

# environment.yml contains dependencies that do not conflict with
# python 2.7.  The second environment contains macs2, which requires
# python 2.7.
RUN conda env create -f /environments/atac.yml &&\
    conda env create -f /environments/macs2.yml

SHELL ["/bin/bash", "-c"]

RUN echo "source activate atac" > ~/.bashrc
ENV PATH /usr/local/envs/atac/bin:$PATH

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

COPY scripts /scripts

CMD ["snakemake", "--directory", "/output","--snakefile", "/scripts/Snakefile","--jobs","4"]
