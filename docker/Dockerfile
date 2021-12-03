FROM nvidia/cuda:11.1-base-ubuntu18.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

# update and install dependencies
RUN apt update && \
    apt install -y software-properties-common && \
    add-apt-repository -y ppa:deadsnakes/ppa && \
    apt update && \
    DEBIAN_FRONTEND=noninteractive apt install -y python3.6 python3-pip wget vim autoconf gcc g++ uuid build-essential \
        bzip2 zlib1g-dev libbz2-dev libcurl4-gnutls-dev liblzma-dev libncurses5-dev libncursesw5-dev libssl-dev && \
    apt clean && \
    apt purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip3 install --upgrade pip==21.2.4 && \
    pip3 install ont-pyguppy-client-lib==5.0.11 megalodon==2.3.5

WORKDIR /opt
RUN wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_5.0.11_linux64.tar.gz && \
    tar xvf ont-guppy-cpu_5.0.11_linux64.tar.gz && \
    mv ont-guppy-cpu ont-guppy-cpu_5.0.11 && \
    rm ont-guppy-cpu_5.0.11_linux64.tar.gz && \
    wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_5.0.11_linux64.tar.gz && \
    tar xvf ont-guppy_5.0.11_linux64.tar.gz && \
    mv ont-guppy ont-guppy_5.0.11 && \
    rm ont-guppy_5.0.11_linux64.tar.gz
ENV GUPPY_GPU_DIR="/opt/ont-guppy_5.0.11"
ENV GUPPY_CPU_DIR="/opt/ont-guppy-cpu_5.0.11"

WORKDIR /opt
RUN wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2 && \
    tar -xvf htslib-1.13.tar.bz2 && \
    rm -r /opt/htslib-1.13.tar.bz2 && \
    cd htslib-1.13 && \
    ./configure && \
    make install
WORKDIR /opt
RUN wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && \
    tar xvf samtools-1.13.tar.bz2 && \
    rm -r /opt/samtools-1.13.tar.bz2 && \
    cd samtools-1.13/ && \
    autoheader && \
    autoconf -Wno-header && \
    ./configure && \
    make install
WORKDIR /opt
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN wget https://github.com/arq5x/bedtools2/archive/refs/tags/v2.30.0.tar.gz && \
    tar xvf v2.30.0.tar.gz && \
    rm -r v2.30.0.tar.gz && \
    cd bedtools2-2.30.0/ && \
    make install
