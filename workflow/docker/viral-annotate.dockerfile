#### viral-annotate.smk #####

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="00307cc41772bade4ec4d37d65038f7bbf27cb1d60e570a9e0d7437231175443"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/eggnog-mapper.yml
#   prefix: /conda-envs/158386e28c523368fcaf484068a6b153
#   name: eggnog-mapper
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - eggnog-mapper=2.0.1=py_1
RUN mkdir -p /conda-envs/158386e28c523368fcaf484068a6b153
COPY workflow/envs/eggnog-mapper.yml /conda-envs/158386e28c523368fcaf484068a6b153/environment.yaml

# Conda environment:
#   source: workflow/envs/metacerberus.yml
#   prefix: /conda-envs/5b969d8fbdf9ab9c8554319f85c9f1e1
#   name: metacerberus
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - fastp=1.3.3
#     - fastqc=0.12.1
#     - metacerberus=1.4.0
#     - pyhmmer=0.10.15
#     - pyrodigal=3.7.1
#     - pyrodigal-gv=0.3.2
#     - python=3.12.13
#     - samtools=1.23.1
#     - scikit-learn=1.9.0
#     - scipy=1.17.1
#     - python-kaleido=0.2.1
RUN mkdir -p /conda-envs/5b969d8fbdf9ab9c8554319f85c9f1e1
COPY workflow/envs/metacerberus.yml /conda-envs/5b969d8fbdf9ab9c8554319f85c9f1e1/environment.yaml

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
#   source: workflow/envs/pharokka.yml
#   prefix: /conda-envs/70f56d9258420de64388e2b8d1e12709
#   name: pharokka
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - pharokka=1.9.1
#     - python=3.11.15
RUN mkdir -p /conda-envs/70f56d9258420de64388e2b8d1e12709
COPY workflow/envs/pharokka.yml /conda-envs/70f56d9258420de64388e2b8d1e12709/environment.yaml

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

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/158386e28c523368fcaf484068a6b153 --file /conda-envs/158386e28c523368fcaf484068a6b153/environment.yaml && \
    conda env create --prefix /conda-envs/5b969d8fbdf9ab9c8554319f85c9f1e1 --file /conda-envs/5b969d8fbdf9ab9c8554319f85c9f1e1/environment.yaml && \
    conda env create --prefix /conda-envs/f24788767f0fd430ce4277dd62814291 --file /conda-envs/f24788767f0fd430ce4277dd62814291/environment.yaml && \
    conda env create --prefix /conda-envs/70f56d9258420de64388e2b8d1e12709 --file /conda-envs/70f56d9258420de64388e2b8d1e12709/environment.yaml && \
    conda env create --prefix /conda-envs/de8c53e37e45136041fa44ef6361f0bb --file /conda-envs/de8c53e37e45136041fa44ef6361f0bb/environment.yaml && \
    conda clean --all -y
