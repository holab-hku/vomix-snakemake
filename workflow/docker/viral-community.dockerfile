#### viral-community.smk #####

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="9250ff5386a028dc6a9ab8f954ef9b3f4a60b4b3d08100de25c25e8d1f7eaf17"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/coverm.yml
#   prefix: /conda-envs/b956c53fb5bdee56d5de5de286af7be1
#   name: coverm
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - bwa=0.7.17=he4a0461_11
#     - coverm=0.6.1=h07ea13f_6
#     - fastani=1.34=h4dfc31f_1
#     - minimap2=2.26=he4a0461_2
#     - samtools=1.9=h10a08f8_12
RUN mkdir -p /conda-envs/b956c53fb5bdee56d5de5de286af7be1
COPY workflow/envs/coverm.yml /conda-envs/b956c53fb5bdee56d5de5de286af7be1/environment.yaml

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

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/b956c53fb5bdee56d5de5de286af7be1 --file /conda-envs/b956c53fb5bdee56d5de5de286af7be1/environment.yaml && \
    conda env create --prefix /conda-envs/4d4b609147c6542c131eb3d54547f609 --file /conda-envs/4d4b609147c6542c131eb3d54547f609/environment.yaml && \
    conda clean --all -y
