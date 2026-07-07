#### viral-host.smk #####

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="9b26b766a5809f3ba437c09bc41c0546e5715bb41e412ed53f94f2c325ad30f8"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/phabox2.yml
#   prefix: /conda-envs/f24788767f0fd430ce4277dd62814291
#   name: phabox2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - phabox=2.1.5=pyhdfd78af_0
RUN mkdir -p /conda-envs/f24788767f0fd430ce4277dd62814291
COPY workflow/envs/phabox2.yml /conda-envs/f24788767f0fd430ce4277dd62814291/environment.yaml

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

RUN conda env create --prefix /conda-envs/f24788767f0fd430ce4277dd62814291 --file /conda-envs/f24788767f0fd430ce4277dd62814291/environment.yaml && \
    conda env create --prefix /conda-envs/6b8a84643dc3ec104144189c2529c0a7 --file /conda-envs/6b8a84643dc3ec104144189c2529c0a7/environment.yaml && \
    conda clean --all -y
