#### prok-annotate.smk #####

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="08a78fdba7aed9f68641ba3fa3d700b2593261f1e057ac5b57792721e82a63ac"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/biobakery3.yml
#   prefix: /conda-envs/7caab56caf267c5dfd6ebaf33db2bb14
#   name: biobakery3
#   channels:
#     - biobakery
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - humann=3.9=py_0
#     - python=3.7.12=hf930737_100_cpython
#     - metaphlan=4.0.6=pyhca03a8a_0
#     - bowtie2=2.5.4=h7071971_4
#     - diamond=2.1.10=h43eeafb_2
RUN mkdir -p /conda-envs/7caab56caf267c5dfd6ebaf33db2bb14
COPY workflow/envs/biobakery3.yml /conda-envs/7caab56caf267c5dfd6ebaf33db2bb14/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/7caab56caf267c5dfd6ebaf33db2bb14 --file /conda-envs/7caab56caf267c5dfd6ebaf33db2bb14/environment.yaml && \
    conda clean --all -y
