FROM ubuntu:16.04
LABEL author="onur.yukselen@umassmed.edu"  description="Docker image containing all requirements for the tRAX pipeline"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH


RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git libtbb-dev g++
    
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc
    
# configure image 
RUN apt-get -y update 
#RUN apt-get -y install software-properties-common build-essential bc
#RUN apt-get -y update

RUN apt-get update && apt-get install -y unzip 
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.0.30.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install
RUN aws --version

RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN conda update -n base -c defaults conda
# Set working directory and copy TRAX software into docker container
COPY . /opt/trax/
ENV PATH /opt/trax:$PATH
RUN conda env create -f /opt/trax/trax_env.yaml
RUN mkdir /rnadb && \
     chmod -R 777 /rnadb && \
     chmod -R 777 /home

RUN conda clean -a
RUN mkdir -p /project /nl /mnt /share
ENV PATH /opt/conda/envs/trax_env/bin:$PATH

