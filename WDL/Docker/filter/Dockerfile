#############################################################
# Dockerfile to build a basic linux environment for bioinformatics
#############################################################

# Set the base image to Ubuntu
FROM ubuntu:20.04

################## METADATA ######################
LABEL base.image="ubuntu:20.04"
LABEL version="1"
LABEL software="enHUNTER base image"
LABEL about.summary="Base image for enHUNTER pipeline"
LABEL about.home="http://bioinformatics.bwh.harvard.edu"
LABEL about.tags="Bioinformatics,Genomics,Transcriptomics,General"
LABEL maintainer="rwang29@bwh.harvard.edu"

# Setup packages
USER root
ENV DEBIAN_FRONTEND noninteractive

RUN mv /etc/apt/sources.list /etc/apt/sources.list.bkp && \
    bash -c 'echo -e "deb mirror://mirrors.ubuntu.com/mirrors.txt xenial main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-updates main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-backports main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-security main restricted universe multiverse\n\n" > /etc/apt/sources.list' && \
    cat /etc/apt/sources.list.bkp >> /etc/apt/sources.list && \
    cat /etc/apt/sources.list

#cleans out the cache 
RUN apt-get clean all && \   
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends apt-utils \
        sudo            \
        bc              \
        autotools-dev   \
        automake        \
        cmake           \
        curl            \
        grep            \
        sed             \
        dpkg            \
        fuse            \
        git             \
        wget            \
        zip             \
        openjdk-8-jre   \
        build-essential \
        pkg-config      \
        python3.6       \
	    python3-dev      \
        python3-pip      \
        bzip2           \
        ca-certificates \
        libglib2.0-0    \
        libssl1.0.0     \
        libssl-dev      \
        libxext6        \
        libsm6          \
        libxrender1     \
        zlib1g-dev &&   \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.0.5-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
 
RUN TINI_VERSION=`curl -L https://github.com/krallin/tini/releases/latest | grep -o "v[0-9]*\.[0-9]*\.[0-9]*" | tail -n1 ` && \
    curl -L "https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini_${TINI_VERSION#?}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

RUN mkdir /data /config

# Add user biodocker with password biodocker
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse biodocker && \
    echo `echo "biodocker\nbiodocker\n" | passwd biodocker` && \
    chown biodocker:biodocker /data && \
    chown biodocker:biodocker /config

# give write permissions to conda folder
RUN chmod 777 -R /opt/conda/

# Change user
USER biodocker

ENV PATH=$PATH:/opt/conda/bin
ENV PATH=$PATH:/home/biodocker/bin
ENV HOME=/home/biodocker

RUN mkdir /home/biodocker/bin

RUN conda config --add channels r
RUN conda config --add channels bioconda

RUN conda upgrade conda

# install R and R packages
RUN conda install -c r r

RUN R -e "install.packages('fitdistrplus',dependencies=TRUE, repos='http://cran.rstudio.com/')" 

# install bedtools
RUN conda install -c bioconda bedtools
RUN conda install -c bioconda samtools
RUN conda install -c bioconda htslib
RUN conda install -c bioconda ucsc-bedgraphtobigwig ucsc-rowstocols ucsc-bigwigmerge

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD ["/bin/bash"]
