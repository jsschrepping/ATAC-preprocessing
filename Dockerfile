# Set the base image to debian based miniconda3 
FROM conda/miniconda3

# File Author / Maintainer 
MAINTAINER JAY-doq

RUN apt-get clean &&\ 
	apt-get update &&\
	apt-get install -y locales

RUN echo "LC_ALL=en_US.UTF-8" >> /etc/environment 
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen 
RUN echo "LANG=en_US.UTF-8" > /etc/locale.conf 
RUN locale-gen en_US.UTF-8

RUN conda update --y conda &&\ 
	conda config --add channels r &&\
	conda config --add channels conda-forge &&\ 
	conda config --add channels bioconda &&\
	conda install -y python=2.7 trimmomatic=0.36 bowtie2=2.3.4.1 picard=2.18.4\
	bedtools=2.27.1 homer=4.9.1 samtools=1.8 macs2=2.1.1.20160309