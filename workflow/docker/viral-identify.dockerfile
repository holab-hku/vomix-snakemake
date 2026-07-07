#### viral-identify.smk #####

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="d34fedf55f8939ca4811245fa6cf3b25846a0888571f75727cd9731031f2283d"

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
#   source: workflow/envs/genomad.yml
#   prefix: /conda-envs/7173b31b5691e611d5cfe5bfa202934b
#   name: genomad
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - genomad=1.8.1=pyhdfd78af_0
#     - keras=2.13.1=pyhd8ed1ab_0
#     - tensorflow=2.13.1=cpu_py310hd1aba9c_1
RUN mkdir -p /conda-envs/7173b31b5691e611d5cfe5bfa202934b
COPY workflow/envs/genomad.yml /conda-envs/7173b31b5691e611d5cfe5bfa202934b/environment.yaml

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

# Conda environment:
#   source: workflow/envs/seqkit-biopython.yml
#   prefix: /conda-envs/6b8a84643dc3ec104144189c2529c0a7
#   name: seqkit
#   channels:
#     - anaconda
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - biopython=1.78=py312h5eee18b_0
#     - numpy=1.26.0=py312hc5e2394_0
#     - pandas=2.1.1=py312h526ad5a_0
#     - python=3.12.0=h996f2a0_0
#     - seqkit=2.6.1=h9ee0642_0
RUN mkdir -p /conda-envs/6b8a84643dc3ec104144189c2529c0a7
COPY workflow/envs/seqkit-biopython.yml /conda-envs/6b8a84643dc3ec104144189c2529c0a7/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/13d30f68c84b09b7da2d01cd98ac951d --file /conda-envs/13d30f68c84b09b7da2d01cd98ac951d/environment.yaml && \
    conda env create --prefix /conda-envs/7173b31b5691e611d5cfe5bfa202934b --file /conda-envs/7173b31b5691e611d5cfe5bfa202934b/environment.yaml && \
    conda env create --prefix /conda-envs/de8c53e37e45136041fa44ef6361f0bb --file /conda-envs/de8c53e37e45136041fa44ef6361f0bb/environment.yaml && \
    conda env create --prefix /conda-envs/57efb9872fbf512b99125d45fc29227f --file /conda-envs/57efb9872fbf512b99125d45fc29227f/environment.yaml && \
    conda env create --prefix /conda-envs/6b8a84643dc3ec104144189c2529c0a7 --file /conda-envs/6b8a84643dc3ec104144189c2529c0a7/environment.yaml && \
    conda clean --all -y
