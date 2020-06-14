# =====================================================================
#  Copied from: https://github.com/erikwaskiewicz/bamsurgeon-nextflow
#  erikwaskiewicz/bamsurgeon Dockerfile
# =====================================================================

# Use miniconda parent image - debian based
FROM continuumio/miniconda:4.7.10
LABEL maintainer="Daniel E. Cook"

# Set the working directory & copy across current directory
# Set python so that it doesn't create .pyc files
WORKDIR /var/app
COPY . /var/app
ENV PYTHONDONTWRITEBYTECODE=true

# Install awk and ps for use with Nextflow reports
# remove apt cache after running to reduce image size
RUN apt-get update \
 && apt-get install --yes --no-install-recommends \
        build-essential \
        gawk \
        git \
        procps \
 && rm -rf /var/lib/apt/lists/*

# Update conda, make the conda environment, clean unnecessary files
# afterwards
RUN conda install --yes --freeze-installed \
        python=3.6 \
        bioconda::bcftools \
        bioconda::bwa \
        bioconda::exonerate \
        bioconda::htslib \
        bioconda::picard=2.22.3-0 \
        bioconda::pysam \
        bioconda::samtools \
        bioconda::velvet \
        conda-forge::numpy \
        conda-forge::scipy \
 && conda clean -afy \
 && find /opt/conda/ -follow -type f -name '*.a' -delete \
 && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
 && find /opt/conda/ -follow -type f -name '*.js.map' -delete

# Download and setup BAMSurgeon - change to last commit before it was
# switched to Python3
RUN git clone https://github.com/adamewing/bamsurgeon.git \
    && cd bamsurgeon \
    && git checkout 20d431ebf5b0101ea2327b84cd520e3faebb4be4 \
    && python setup.py install

