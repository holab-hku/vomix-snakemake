#### prok-community.smk #####

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="25687c7e55be4b65411804a596ca205e0ab70b20811a1d4e92b4f331fc5b88f9"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/metaphlan.yml
#   prefix: /conda-envs/cf196b66330ad8bdfa1e58ede220c9b0
#   name: metaphlan
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - metaphlan=4.0.6=pyhca03a8a_0
#     - python=3.10.0
RUN mkdir -p /conda-envs/cf196b66330ad8bdfa1e58ede220c9b0
COPY workflow/envs/metaphlan.yml /conda-envs/cf196b66330ad8bdfa1e58ede220c9b0/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/cf196b66330ad8bdfa1e58ede220c9b0 --file /conda-envs/cf196b66330ad8bdfa1e58ede220c9b0/environment.yaml && \
    conda clean --all -y
