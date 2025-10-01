# Adapted from https://github.com/ding-lab/HotspotFilter
FROM continuumio/miniconda3

MAINTAINER Matthew Wyczalkowski <m.wyczalkowski@wustl.edu>

RUN conda create -n env python=3.6
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH

WORKDIR /opt
# https://github.com/geerlingguy/ansible-role-java/issues/64
RUN mkdir -p /usr/share/man/man1

# run update
#RUN apt-get update -y && apt-get install -y \
#    curl \
#    libnss-sss \
#    vim

RUN conda install -c bioconda scipy bioconda::pysam conda-forge::numpy bioconda::pyvcf3 conda-forge::scikit-learn bioconda::cigar 

RUN conda install -c bioconda cutesv

# https://github.com/ding-lab/BreakPointSurveyor
#RUN git clone --recursive https://github.com/ding-lab/BreakPointSurveyor.git

#COPY . /opt/HotspotFilter
