FROM ubuntu:xenial

LABEL author="Matt Wyczalkowski" \
      maintainer="m.wyczalkowski@wustl.edu"

RUN apt-get update && apt-get install -y \
    build-essential \
    cpanminus \
    locales \
    bzip2 \
    tabix \
    && apt-get clean

# Set the locale so perl doesn't complain
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

CMD ["/bin/bash"]
