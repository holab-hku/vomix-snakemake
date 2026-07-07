#### viral-host.smk (—iphop-host) #####

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="736022dafde9d7f12683c9254c1a20f358737710f9a993d863394ffd054f7f9f"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/iphop.yml
#   prefix: /conda-envs/bac0f064995993010230def5ad9686bf
#   name: iphop
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - iphop=1.3.3=pyhdfd78af_0
#     - python=3.8.15=h257c98d_0_cpython
RUN mkdir -p /conda-envs/bac0f064995993010230def5ad9686bf
COPY workflow/envs/iphop.yml /conda-envs/bac0f064995993010230def5ad9686bf/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/bac0f064995993010230def5ad9686bf --file /conda-envs/bac0f064995993010230def5ad9686bf/environment.yaml && \
    conda clean --all -y
