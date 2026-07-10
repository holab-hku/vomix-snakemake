#### viral-annotate.smk #####

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="7e45ffa61eabd33ca97daeca935a05a96fde64de06a1ce82c4e0fc1cb7e071da"
# Conda environment:
#   source: workflow/envs/biobakery3.yml
#   prefix: /conda-envs/6c71dbf8184f5937b984da1448b7a1cc
#   name: biobakery3
#   channels:
#     - biobakery
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - humann=3.9
#     - python=3.7.12
#     - metaphlan=4.0.6
#     - bowtie2=2.5.4
#     - diamond=2.1.10
RUN mkdir -p /conda-envs/6c71dbf8184f5937b984da1448b7a1cc
COPY workflow/envs/biobakery3.yml /conda-envs/6c71dbf8184f5937b984da1448b7a1cc/environment.yaml
# Conda environment:
#   source: workflow/envs/checkm2.yml
#   prefix: /conda-envs/8cd1f42070afc8f3535013e390ce4e5b
#   name: checkm2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - checkm2=1.0.2
#     - prodigal=2.6.3
#     - scikit-learn=0.23.2
#     - scipy=1.8.0
#     - tensorflow=2.4.0
#     - tensorflow-base=2.4.0
RUN mkdir -p /conda-envs/8cd1f42070afc8f3535013e390ce4e5b
COPY workflow/envs/checkm2.yml /conda-envs/8cd1f42070afc8f3535013e390ce4e5b/environment.yaml
# Conda environment:
#   source: workflow/envs/checkv.yml
#   prefix: /conda-envs/c05eb9865829d1fa8c0fe1d9af4d4bea
#   name: checkv
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - checkv=1.0.1
#     - blast=2.15.0
#     - diamond=2.1.8
#     - hmmer=3.4
#     - python=3.10.13
RUN mkdir -p /conda-envs/c05eb9865829d1fa8c0fe1d9af4d4bea
COPY workflow/envs/checkv.yml /conda-envs/c05eb9865829d1fa8c0fe1d9af4d4bea/environment.yaml
# Conda environment:
#   source: workflow/envs/dram.yml
#   prefix: /conda-envs/bfbcf982e543752be0db295d2d8fc7ca
#   name: dram
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python=3.10
#     - pandas=1.5.2
#     - pytest=7.2.0
#     - scikit-bio=0.5.7
#     - prodigal=2.6.3
#     - mmseqs2==13.45111
#     - hmmer=3.3.2
#     - trnascan-se=2.0.11
#     - scipy=1.8.1
#     - sqlalchemy=1.4.46
#     - barrnap=0.9
#     - altair=4.2.0
#     - openpyxl=3.0.10
#     - networkx=2.8.8
#     - ruby=3.1.2
#     - parallel=20221122
#     - pip
#     - pip:
#       - DRAM-bio==1.5.0
RUN mkdir -p /conda-envs/bfbcf982e543752be0db295d2d8fc7ca
COPY workflow/envs/dram.yml /conda-envs/bfbcf982e543752be0db295d2d8fc7ca/environment.yaml
# Conda environment:
#   source: workflow/envs/eggnog-mapper.yml
#   prefix: /conda-envs/b8be70dcd988c91ea813792b098e3fd9
#   name: eggnog-mapper
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - eggnog-mapper=2.0.1
RUN mkdir -p /conda-envs/b8be70dcd988c91ea813792b098e3fd9
COPY workflow/envs/eggnog-mapper.yml /conda-envs/b8be70dcd988c91ea813792b098e3fd9/environment.yaml
# Conda environment:
#   source: workflow/envs/genomad.yml
#   prefix: /conda-envs/303cc820d42025793664c0fcc95e1147
#   name: genomad
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - genomad=1.8.1
#     - keras=2.13.1
#     - tensorflow=2.13.1
RUN mkdir -p /conda-envs/303cc820d42025793664c0fcc95e1147
COPY workflow/envs/genomad.yml /conda-envs/303cc820d42025793664c0fcc95e1147/environment.yaml
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
#   source: workflow/envs/iphop.yml
#   prefix: /conda-envs/a4751a6f51373a21124424968cc15c53
#   name: iphop
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - iphop=1.3.3
#     - python=3.8.15
RUN mkdir -p /conda-envs/a4751a6f51373a21124424968cc15c53
COPY workflow/envs/iphop.yml /conda-envs/a4751a6f51373a21124424968cc15c53/environment.yaml
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
#   prefix: /conda-envs/d4d369c29c6d48bdb1c3ed652ca3aed1
#   name: phabox2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - phabox=2.1.5
RUN mkdir -p /conda-envs/d4d369c29c6d48bdb1c3ed652ca3aed1
COPY workflow/envs/phabox2.yml /conda-envs/d4d369c29c6d48bdb1c3ed652ca3aed1/environment.yaml
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
#   prefix: /conda-envs/3718f9432b895f1c61eb898c6c547322
#   name: prodigal-gv
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - prodigal-gv=2.11.0
RUN mkdir -p /conda-envs/3718f9432b895f1c61eb898c6c547322
COPY workflow/envs/prodigal-gv.yml /conda-envs/3718f9432b895f1c61eb898c6c547322/environment.yaml
# Conda environment:
#   source: workflow/envs/vibrant.yml
#   prefix: /conda-envs/562ec28b54b50c320491196954064d79
#   name: vibrant
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - biopython=1.79
#       - hmmer=3.4
#       - matplotlib=3.3.2
#       - numpy=1.19.5
#       - pandas=0.25.3
#       - prodigal=2.6.3
#       - python=3.7.12
#       - scikit-learn=0.21.3
#       - scipy=1.7.3
#       - seaborn=0.12.2
#       - tk=8.6.13
#       - vibrant=1.2.1
RUN mkdir -p /conda-envs/562ec28b54b50c320491196954064d79
COPY workflow/envs/vibrant.yml /conda-envs/562ec28b54b50c320491196954064d79/environment.yaml
# Conda environment:
#   source: workflow/envs/virsorter2.yml
#   prefix: /conda-envs/368c0dfb9ae7365d73597249b1529bba
#   name: virsorter2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - virsorter=2.2.4
RUN mkdir -p /conda-envs/368c0dfb9ae7365d73597249b1529bba
COPY workflow/envs/virsorter2.yml /conda-envs/368c0dfb9ae7365d73597249b1529bba/environment.yaml

RUN conda env create --prefix /conda-envs/6c71dbf8184f5937b984da1448b7a1cc --file /conda-envs/6c71dbf8184f5937b984da1448b7a1cc/environment.yaml && \
    conda env create --prefix /conda-envs/8cd1f42070afc8f3535013e390ce4e5b --file /conda-envs/8cd1f42070afc8f3535013e390ce4e5b/environment.yaml && \
    conda env create --prefix /conda-envs/c05eb9865829d1fa8c0fe1d9af4d4bea --file /conda-envs/c05eb9865829d1fa8c0fe1d9af4d4bea/environment.yaml && \
    conda env create --prefix /conda-envs/bfbcf982e543752be0db295d2d8fc7ca --file /conda-envs/bfbcf982e543752be0db295d2d8fc7ca/environment.yaml && \
    conda env create --prefix /conda-envs/b8be70dcd988c91ea813792b098e3fd9 --file /conda-envs/b8be70dcd988c91ea813792b098e3fd9/environment.yaml && \
    conda env create --prefix /conda-envs/303cc820d42025793664c0fcc95e1147 --file /conda-envs/303cc820d42025793664c0fcc95e1147/environment.yaml && \
    conda env create --prefix /conda-envs/3094ed1ebef010ad0fae5d98c40a1a95 --file /conda-envs/3094ed1ebef010ad0fae5d98c40a1a95/environment.yaml && \
    conda env create --prefix /conda-envs/4ecab50c611bd20ca73e77590bda0e41 --file /conda-envs/4ecab50c611bd20ca73e77590bda0e41/environment.yaml && \
    conda env create --prefix /conda-envs/a4751a6f51373a21124424968cc15c53 --file /conda-envs/a4751a6f51373a21124424968cc15c53/environment.yaml && \
    conda env create --prefix /conda-envs/5b969d8fbdf9ab9c8554319f85c9f1e1 --file /conda-envs/5b969d8fbdf9ab9c8554319f85c9f1e1/environment.yaml && \
    conda env create --prefix /conda-envs/d4d369c29c6d48bdb1c3ed652ca3aed1 --file /conda-envs/d4d369c29c6d48bdb1c3ed652ca3aed1/environment.yaml && \
    conda env create --prefix /conda-envs/70f56d9258420de64388e2b8d1e12709 --file /conda-envs/70f56d9258420de64388e2b8d1e12709/environment.yaml && \
    conda env create --prefix /conda-envs/3718f9432b895f1c61eb898c6c547322 --file /conda-envs/3718f9432b895f1c61eb898c6c547322/environment.yaml && \
    conda env create --prefix /conda-envs/562ec28b54b50c320491196954064d79 --file /conda-envs/562ec28b54b50c320491196954064d79/environment.yaml && \
    conda env create --prefix /conda-envs/368c0dfb9ae7365d73597249b1529bba --file /conda-envs/368c0dfb9ae7365d73597249b1529bba/environment.yaml && \
    conda clean --all -y
