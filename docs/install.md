# Installation

## Conda & Mamba

You can install vomix-snakemake in you computer using general-purpose package managers ( `mamba`, `conda`). .

::::{tab-set}

:::{tab-item} Conda
[Conda](https://docs.conda.io/projects/conda/en/stable/) is a package manager that handles all your dependencies for you. To install vomix-snakemake using Conda, you can create the environment from the repository environment file.

```bash
# Download GitHub directory
git clone https://github.com/holab-hku/vomix-snakemake.git
cd vomix-snakemake

# Install base environment
conda env create -f environment.yml

# Activate environment
conda activate vomix

# Verify Installation
snakemake -v
```

:::

:::{tab-item} Mamba
[Mamba](https://mamba.readthedocs.io/en/latest/) is a package manager that handles all your dependencies for you. To install vomix-snakemake using Mamba, you can create the environment directly from the repository environment file.

```bash
# Download GitHub directory
git clone https://github.com/holab-hku/vomix-snakemake.git
cd vomix-snakemake

# Install base environment
mamba env create -f environment.yml

# Activate environment
mamba activate vomix

# Verify Installation
snakemake -v
```

:::

:::{tab-item} Conda Lock
[Conda Lock](https://github.com/conda/conda-lock) provides a pinned dependency lock file that is often the most reliable fallback when other installation methods fail.

```bash
# Download GitHub directory
git clone https://github.com/holab-hku/vomix-snakemake.git
cd vomix-snakemake

# Create a conda-lock environment 
# if you don't want to install into base environment
conda create -n conda-lock -c conda-forge conda-lock=4.0.2 -y
conda activate conda-lock

# Use conda-lock to install the environment from the repository lock file
conda-lock install --name vomix conda-lock.yml
conda deactivate # deactive conda-lock environment
conda activate vomix # activate vomix environment

# Verify Installation
snakemake -v
```

:::

::::

```{admonition} Conda update
:class: attention
If you are using conda or mamba, make sure your conda installation is up to date before proceeding. You can update it with `conda update -n base -c defaults conda`
```

```{admonition} Conda Channel Priorities
:class: attention
If you are using conda or mamba, make sure to set channel orders correctly and set channel priority to strict. Via the `conda config --add channels defaults`, `conda config --add channels bioconda`, `conda config --add channels conda-forge`, and `conda config --set channel_priority strict` respectively. For mamba replace `conda` with `mamba` respectively.  
```

```{admonition} Conda Lock as a fallback
:class: note
If the standard conda or mamba installation methods do not work, `conda-lock` is usually the best approach because it installs the pinned dependency set from the repository lock file. Snakemake and many bioinformatics tools rely heavily on POSIX-compliant workflows, so they are currently best supported on non-Windows operating systems.
```

## Docker, Apptainer, and Singularity

vomix-snakemake is built on a snakemake back-end, which facilitates native `Docker` deployment via a `Apptainer` (formerly `Singularity`) `.sif` image. The container image generated contains explicitly each conda environment mounted on top of a base operating system. Containers are preferred for the most robus forms of reproducibility, whereas `conda` and `mamba` installations might not work on Windows or Mac-ARM systems.

::::{tab-set}
:::{tab-item} Apptainer/Singularity

```bash
snakemake --software-deployment-method apptainer
```

:::

::::

```{admonition} Using Conda within Containers
:class: note
Each rule in vomix-snakemake's underlying snakemake files depends on a speicifc conda environment. To run what Snakemake calls `Ad-hoc combination of Conda package management with containers`, which is essentially running apptainer containers with conda enviornments installed within them, you need to use the `--sdm conda apptainer` option. This allows true full reproducibility. 
```

If you only use `--sdm apptainer` Snakemake will not launch any conda environments and hence jobs will fail. If you use `--sdm apptainer` it will try and re-install conda enviornments in your local `.snakemake/conda` folder, which counteracts the purpose of containers.

```

## {octicon}`book;0.85em` Troubleshooting Guide

We have specific guidelines for troubleshooting vomix-snakemake so we can help you out in your analysis journey as efficiently as possible! If you run into any unexpected errors, warnings, etc. please visit our [Troubleshooting Guide](/troubleshoot.md).

## {octicon}`bug;0.85em` Report a bug to us

Have any questions or you've found a bug during your analysis? Please don't hesitate to report it to us by making an issue on our [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/holab-hku/vOMIX-MEGA/issues/new).
