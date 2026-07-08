#### assembly.smk #####

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="7e8c09c67a6d243bd5e0719a666e08e325fb0a067f8659efabd8c6c4a792a56e"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/megahit.yml
#   prefix: /conda-envs/cb1b9291fea5d2dd0747641376bbb1f8
#   name: megahit
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - megahit=1.2.9=h43eeafb_4
RUN mkdir -p /conda-envs/cb1b9291fea5d2dd0747641376bbb1f8
COPY workflow/envs/megahit.yml /conda-envs/cb1b9291fea5d2dd0747641376bbb1f8/environment.yaml

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

RUN conda env create --prefix /conda-envs/cb1b9291fea5d2dd0747641376bbb1f8 --file /conda-envs/cb1b9291fea5d2dd0747641376bbb1f8/environment.yaml && \
    conda env create --prefix /conda-envs/6b8a84643dc3ec104144189c2529c0a7 --file /conda-envs/6b8a84643dc3ec104144189c2529c0a7/environment.yaml && \
    conda clean --all -y
