#### preprocess.smk #####

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="f9b26b83619acaa7348e7d623a8a7de1a4224a2712eac7932fb5ba8e33ca3733"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/fastp.yml
#   prefix: /conda-envs/111ef10e87d672fe57ebd1dd2cbe89de
#   name: fastp
#   channels:
#     - anaconda
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - fastp=0.22.0=h2e03b76_0
RUN mkdir -p /conda-envs/111ef10e87d672fe57ebd1dd2cbe89de
COPY workflow/envs/fastp.yml /conda-envs/111ef10e87d672fe57ebd1dd2cbe89de/environment.yaml

# Conda environment:
#   source: workflow/envs/hostile.yml
#   prefix: /conda-envs/4ecab50c611bd20ca73e77590bda0e41
#   name: hostile
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - hostile=2.0.1
RUN mkdir -p /conda-envs/4ecab50c611bd20ca73e77590bda0e41
COPY workflow/envs/hostile.yml /conda-envs/4ecab50c611bd20ca73e77590bda0e41/environment.yaml

# Conda environment:
#   source: workflow/envs/multiqc.yml
#   prefix: /conda-envs/c739858ae050f21e5cf8a089d1aad71e
#   name: multiqc
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - multiqc=1.17=pyhdfd78af_1
RUN mkdir -p /conda-envs/c739858ae050f21e5cf8a089d1aad71e
COPY workflow/envs/multiqc.yml /conda-envs/c739858ae050f21e5cf8a089d1aad71e/environment.yaml

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

# Conda environment:
#   source: workflow/envs/sratools-pigz.yml
#   prefix: /conda-envs/6651a7798ae83189ba4a13938616c001
#   name: sra-tools
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - pigz=2.8=h2797004_0
#     - sra-tools=3.2.1
RUN mkdir -p /conda-envs/6651a7798ae83189ba4a13938616c001
COPY workflow/envs/sratools-pigz.yml /conda-envs/6651a7798ae83189ba4a13938616c001/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/111ef10e87d672fe57ebd1dd2cbe89de --file /conda-envs/111ef10e87d672fe57ebd1dd2cbe89de/environment.yaml && \
    conda env create --prefix /conda-envs/4ecab50c611bd20ca73e77590bda0e41 --file /conda-envs/4ecab50c611bd20ca73e77590bda0e41/environment.yaml && \
    conda env create --prefix /conda-envs/c739858ae050f21e5cf8a089d1aad71e --file /conda-envs/c739858ae050f21e5cf8a089d1aad71e/environment.yaml && \
    conda env create --prefix /conda-envs/6b8a84643dc3ec104144189c2529c0a7 --file /conda-envs/6b8a84643dc3ec104144189c2529c0a7/environment.yaml && \
    conda env create --prefix /conda-envs/6651a7798ae83189ba4a13938616c001 --file /conda-envs/6651a7798ae83189ba4a13938616c001/environment.yaml && \
    conda clean --all -y
