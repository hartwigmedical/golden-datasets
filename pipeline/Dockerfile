# Set the base image to
FROM google/cloud-sdk:latest
ADD Miniconda3-latest-Linux-x86_64.sh Miniconda3-latest-Linux-x86_64.sh
ADD hmf-nextflow-docker hmf-nextflow-docker
RUN yes yes | bash Miniconda3-latest-Linux-x86_64.sh
RUN curl -s https://get.nextflow.io | bash

ENTRYPOINT ./hmf-nextflow-docker
