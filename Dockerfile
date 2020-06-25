FROM python:3.8-slim

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
    build-essential \
    cmake \
    gfortran \
    libblas-dev \
    libbz2-dev \
    liblapack-dev \
    make \
    unzip \
    wget \
    zlib1g-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /app

## Set SVDB version
ARG SVDB_VERSION=2.4.0

## Add some info
LABEL base_image="python:3.8-slim"
LABEL software="svdb"
LABEL software.version=${SVDB_VERSION}

RUN pip install --no-cache-dir cython numpy

## Download and install
RUN wget https://github.com/J35P312/SVDB/archive/${SVDB_VERSION}.zip && \
    unzip ${SVDB_VERSION}.zip && \
    cd SVDB-${SVDB_VERSION} && \
    pip install -e . && \
    rm /app/${SVDB_VERSION}.zip

ENTRYPOINT ["svdb"]
CMD ["--help"]
