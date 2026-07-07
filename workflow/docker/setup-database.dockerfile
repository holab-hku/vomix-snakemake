#### setup-database.smk #####

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="096de052ae56f327a7c4828f89db25827413e3f70a92544e113b13eab64805d1"

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
#   source: workflow/envs/gtdbtk.yml
#   prefix: /conda-envs/3094ed1ebef010ad0fae5d98c40a1a95
#   name: gtdbtk
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - gtdbtk=2.7.2
RUN mkdir -p /conda-envs/3094ed1ebef010ad0fae5d98c40a1a95
COPY workflow/envs/gtdbtk.yml /conda-envs/3094ed1ebef010ad0fae5d98c40a1a95/environment.yaml

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

RUN conda env create --prefix /conda-envs/13d30f68c84b09b7da2d01cd98ac951d --file /conda-envs/13d30f68c84b09b7da2d01cd98ac951d/environment.yaml && \
    conda env create --prefix /conda-envs/7173b31b5691e611d5cfe5bfa202934b --file /conda-envs/7173b31b5691e611d5cfe5bfa202934b/environment.yaml && \
    conda env create --prefix /conda-envs/3094ed1ebef010ad0fae5d98c40a1a95 --file /conda-envs/3094ed1ebef010ad0fae5d98c40a1a95/environment.yaml && \
    conda env create --prefix /conda-envs/f24788767f0fd430ce4277dd62814291 --file /conda-envs/f24788767f0fd430ce4277dd62814291/environment.yaml && \
    conda clean --all -y
