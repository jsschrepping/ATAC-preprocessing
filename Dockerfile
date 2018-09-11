# Set the base image to debian based miniconda2
FROM conda/miniconda2

# File Author / Maintainer
MAINTAINER JAY-DoQ

COPY environment.yml .

RUN conda env create -f /environment.yml &&\
    source activate atac
