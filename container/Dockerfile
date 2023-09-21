FROM  nvidia/cuda:11.0.3-cudnn8-runtime-ubuntu20.04

# Install update apt and install wget package
RUN apt-get update && \
    apt-get install -y wget procps

# Create user analysis
RUN useradd -m spatial_annotation && \
    mkdir /home/spatial_annotation/.ssh && \
    chown -R spatial_annotation:spatial_annotation /home/spatial_annotation/.ssh

# switch user to analysis
USER spatial_annotation

# Install miniconda
RUN cd /home/spatial_annotation && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /home/spatial_annotation/miniconda3 && \
    rm Miniconda3-latest-Linux-x86_64.sh

# Add miniconda to path
ENV PATH=/home/spatial_annotation/miniconda3/bin:$PATH


# Update base conda env with conda-env yaml file
COPY conda.yml /home/spatial_annotation/environment.yml
RUN conda env update -n base -f /home/spatial_annotation/environment.yml