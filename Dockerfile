# Set the base image to debian based miniconda3 
FROM conda/miniconda3

# File Author / Maintainer 
MAINTAINER JAY-doq

RUN apt-get clean &&\ 
    apt-get update &&\
    apt-get install -y locales

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

RUN conda update -y conda

COPY scripts/environment.yml environment.yml

RUN conda env create -f /environment.yml &&\
    source activate atac
