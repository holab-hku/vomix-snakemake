#### checkv-pyhmmer.smk #####

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="1018ad4a93af5008c91d52c86173b6075e4be1d310a095e1c9ce358b410891c1"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/checkv.yml
#   prefix: /conda-envs/13d30f68c84b09b7da2d01cd98ac951d
#   name: checkv
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - checkv=1.0.1=pyhdfd78af_0
#     - blast=2.15.0
#     - diamond=2.1.8
#     - hmmer=3.4
#     - python=3.10.13
RUN mkdir -p /conda-envs/13d30f68c84b09b7da2d01cd98ac951d
COPY workflow/envs/checkv.yml /conda-envs/13d30f68c84b09b7da2d01cd98ac951d/environment.yaml

# Conda environment:
#   source: workflow/envs/prodigal-gv.yml
#   prefix: /conda-envs/de8c53e37e45136041fa44ef6361f0bb
#   name: prodigal-gv
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - prodigal-gv=2.11.0=he4a0461_2
RUN mkdir -p /conda-envs/de8c53e37e45136041fa44ef6361f0bb
COPY workflow/envs/prodigal-gv.yml /conda-envs/de8c53e37e45136041fa44ef6361f0bb/environment.yaml

# Conda environment:
#   source: workflow/envs/pyhmmer.yml
#   prefix: /conda-envs/57efb9872fbf512b99125d45fc29227f
#   name: pyhmmer
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - pyhmmer=0.10.4=py310h4b81fae_0
RUN mkdir -p /conda-envs/57efb9872fbf512b99125d45fc29227f
COPY workflow/envs/pyhmmer.yml /conda-envs/57efb9872fbf512b99125d45fc29227f/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/13d30f68c84b09b7da2d01cd98ac951d --file /conda-envs/13d30f68c84b09b7da2d01cd98ac951d/environment.yaml && \
    conda env create --prefix /conda-envs/de8c53e37e45136041fa44ef6361f0bb --file /conda-envs/de8c53e37e45136041fa44ef6361f0bb/environment.yaml && \
    conda env create --prefix /conda-envs/57efb9872fbf512b99125d45fc29227f --file /conda-envs/57efb9872fbf512b99125d45fc29227f/environment.yaml && \
    conda clean --all -y
