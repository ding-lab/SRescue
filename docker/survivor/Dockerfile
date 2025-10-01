FROM gcc:4.9

MAINTAINER Matthew Wyczalkowski <m.wyczalkowski@wustl.edu>

WORKDIR /opt

# run update
#RUN apt-get update -y && apt-get install -y \
#    curl \
#    libnss-sss \
#    vim

RUN git clone https://github.com/fritzsedlazeck/SURVIVOR.git && cd SURVIVOR/Debug && make

