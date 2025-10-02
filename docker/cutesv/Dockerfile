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

#RUN conda install -c bioconda scipy bioconda::pysam conda-forge::numpy bioconda::pyvcf3 conda-forge::scikit-learn bioconda::cigar 
RUN conda install -c bioconda scipy pysam numpy pyvcf3 scikit-learn cigar 

#RUN conda install -c bioconda cutesv
#RUN git clone https://github.com/tjiangHIT/cuteSV.git --branch cuteSV-v2.1.1 && cd cuteSV/ && python setup.py install
RUN pip install cuteSV==2.1.1

#RUN git clone https://github.com/tjiangHIT/cuteFC.git && cd cuteFC/ && python setup.py install

# https://github.com/ding-lab/BreakPointSurveyor
#RUN git clone --recursive https://github.com/ding-lab/BreakPointSurveyor.git

#COPY . /opt/HotspotFilter
