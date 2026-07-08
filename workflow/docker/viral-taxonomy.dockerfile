#### viral-taxonomy.smk #####

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="8c5a7732f73b566264df23958369aa1ac8c022b08649c29d4b5955cf272c5e8b"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/ete3.yml
#   prefix: /conda-envs/4d4b609147c6542c131eb3d54547f609
#   name: taxonomy
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - biopython=1.81=py312h98912ed_1
#     - ete3=3.1.3=pyhd8ed1ab_0
#     - numpy=1.26.2=py312heda63a1_0
#     - pandas=2.1.3=py312hfb8ada1_0
#     - python=3.12.0=hab00c5b_0_cpython
#     - scipy=1.11.4=py312heda63a1_0
RUN mkdir -p /conda-envs/4d4b609147c6542c131eb3d54547f609
COPY workflow/envs/ete3.yml /conda-envs/4d4b609147c6542c131eb3d54547f609/environment.yaml

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

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/4d4b609147c6542c131eb3d54547f609 --file /conda-envs/4d4b609147c6542c131eb3d54547f609/environment.yaml && \
    conda env create --prefix /conda-envs/7173b31b5691e611d5cfe5bfa202934b --file /conda-envs/7173b31b5691e611d5cfe5bfa202934b/environment.yaml && \
    conda env create --prefix /conda-envs/f24788767f0fd430ce4277dd62814291 --file /conda-envs/f24788767f0fd430ce4277dd62814291/environment.yaml && \
    conda clean --all -y
