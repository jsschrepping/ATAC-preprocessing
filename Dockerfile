# Set the base image to debian based miniconda2
FROM conda/miniconda3

# File Author / Maintainer
MAINTAINER JAY-DoQ

COPY environments /environments

# environment.yml contains dependencies that do not conflict with
# python 2.7.  The second environment contains macs2, which requires
# python 2.7.
RUN conda env create -f /environments/atac.yml &&\
    conda env create -f /environments/macs2.yml

SHELL ["/bin/bash", "-c"]

RUN source activate atac
