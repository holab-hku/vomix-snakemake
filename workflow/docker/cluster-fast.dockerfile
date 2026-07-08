#### cluster-fast.smk #####

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="fcf9e064e501c17ac3bf960b72ef21b8970a25241fdea959998c6d1b2c29498c"

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
    conda env create --prefix /conda-envs/6b8a84643dc3ec104144189c2529c0a7 --file /conda-envs/6b8a84643dc3ec104144189c2529c0a7/environment.yaml && \
    conda clean --all -y
