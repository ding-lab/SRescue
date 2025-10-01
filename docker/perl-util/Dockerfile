FROM ubuntu:xenial

LABEL author="Matt Wyczalkowski" \
      maintainer="m.wyczalkowski@wustl.edu"

RUN apt-get update && apt-get install -y \
    build-essential \
    cpanminus \
    && apt-get clean

ENV LC_ALL="C"

CMD ["/bin/bash"]
