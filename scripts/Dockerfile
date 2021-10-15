# Use official conda image
FROM conda/miniconda3-centos7

# Set the working directory to /app
WORKDIR /app

# Copy the script directory contents into the container at /app
ADD "https://github.com/EUCANCan/golden-datasets/tree/devel/scripts" /app

# Create conda envs
RUN yum install -y which   \
&& yum clean all \
&& conda env create -n eucancan -f /app/environment_snv.yml \
&& conda env create -n eucancan_sv -f /app/environment.yml \
&& conda env create -n benchmarking_plots -f /app/environment_plots.yml  

ENTRYPOINT [ "bash", "ingest_snv.sh" ]