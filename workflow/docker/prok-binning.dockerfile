#### prok-binning.smk #####

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="f2efa9f235539ddb6d74fbc0cda3fb40f408425223fc764063a97a03c3505a8b"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/checkm2.yml
#   prefix: /conda-envs/be4ed925f0a032f4730b14c0519a7031
#   name: checkm2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - checkm2=1.0.2=pyh7cba7a3_0
#     - prodigal=2.6.3=h7b50bb2_10
#     - scikit-learn=0.23.2=py38h5d63f67_3
#     - scipy=1.8.0=py38h56a6a73_1
#     - tensorflow=2.4.0=py38h578d9bd_0
#     - tensorflow-base=2.4.0=py38h01d9eeb_0
RUN mkdir -p /conda-envs/be4ed925f0a032f4730b14c0519a7031
COPY workflow/envs/checkm2.yml /conda-envs/be4ed925f0a032f4730b14c0519a7031/environment.yaml

# Conda environment:
#   source: workflow/envs/concoct.yml
#   prefix: /conda-envs/c31eace9657ebb7dc76ece21eef79249
#   name: concoct
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - concoct=1.1.0=py312h245ed52_6
#     - python=3.12.6=hc5c86c4_1_cpython
RUN mkdir -p /conda-envs/c31eace9657ebb7dc76ece21eef79249
COPY workflow/envs/concoct.yml /conda-envs/c31eace9657ebb7dc76ece21eef79249/environment.yaml

# Conda environment:
#   source: workflow/envs/dastool.yml
#   prefix: /conda-envs/b30e77ed52893c366d96cbeb6be03b0a
#   name: dastool
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - das_tool=1.1.7=r43hdfd78af_0
RUN mkdir -p /conda-envs/b30e77ed52893c366d96cbeb6be03b0a
COPY workflow/envs/dastool.yml /conda-envs/b30e77ed52893c366d96cbeb6be03b0a/environment.yaml

# Conda environment:
#   source: workflow/envs/galah.yml
#   prefix: /conda-envs/297fee73a3ca64ef4e4a04e526d3f224
#   name: galah
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - galah=0.4.2=h7b50bb2_1
RUN mkdir -p /conda-envs/297fee73a3ca64ef4e4a04e526d3f224
COPY workflow/envs/galah.yml /conda-envs/297fee73a3ca64ef4e4a04e526d3f224/environment.yaml

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
#   source: workflow/envs/maxbin2.yml
#   prefix: /conda-envs/4584f0e5705cb1b9415797d65218aeba
#   name: maxbin2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - maxbin2=2.2.7=hdbdd923_5
RUN mkdir -p /conda-envs/4584f0e5705cb1b9415797d65218aeba
COPY workflow/envs/maxbin2.yml /conda-envs/4584f0e5705cb1b9415797d65218aeba/environment.yaml

# Conda environment:
#   source: workflow/envs/metabat2.yml
#   prefix: /conda-envs/88af96c64a666726623f7325857fdcc6
#   name: metabat2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - metabat2=2.17=hd498684_0
RUN mkdir -p /conda-envs/88af96c64a666726623f7325857fdcc6
COPY workflow/envs/metabat2.yml /conda-envs/88af96c64a666726623f7325857fdcc6/environment.yaml

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
#   source: workflow/envs/strobealign.yml
#   prefix: /conda-envs/12b77e3b8e18f601c6cd28cf83cda940
#   name: strobealign
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - strobealign=0.13.0
#     - samtools=1.23.1
RUN mkdir -p /conda-envs/12b77e3b8e18f601c6cd28cf83cda940
COPY workflow/envs/strobealign.yml /conda-envs/12b77e3b8e18f601c6cd28cf83cda940/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/be4ed925f0a032f4730b14c0519a7031 --file /conda-envs/be4ed925f0a032f4730b14c0519a7031/environment.yaml && \
    conda env create --prefix /conda-envs/c31eace9657ebb7dc76ece21eef79249 --file /conda-envs/c31eace9657ebb7dc76ece21eef79249/environment.yaml && \
    conda env create --prefix /conda-envs/b30e77ed52893c366d96cbeb6be03b0a --file /conda-envs/b30e77ed52893c366d96cbeb6be03b0a/environment.yaml && \
    conda env create --prefix /conda-envs/297fee73a3ca64ef4e4a04e526d3f224 --file /conda-envs/297fee73a3ca64ef4e4a04e526d3f224/environment.yaml && \
    conda env create --prefix /conda-envs/3094ed1ebef010ad0fae5d98c40a1a95 --file /conda-envs/3094ed1ebef010ad0fae5d98c40a1a95/environment.yaml && \
    conda env create --prefix /conda-envs/4584f0e5705cb1b9415797d65218aeba --file /conda-envs/4584f0e5705cb1b9415797d65218aeba/environment.yaml && \
    conda env create --prefix /conda-envs/88af96c64a666726623f7325857fdcc6 --file /conda-envs/88af96c64a666726623f7325857fdcc6/environment.yaml && \
    conda env create --prefix /conda-envs/6b8a84643dc3ec104144189c2529c0a7 --file /conda-envs/6b8a84643dc3ec104144189c2529c0a7/environment.yaml && \
    conda env create --prefix /conda-envs/12b77e3b8e18f601c6cd28cf83cda940 --file /conda-envs/12b77e3b8e18f601c6cd28cf83cda940/environment.yaml && \
    conda clean --all -y
