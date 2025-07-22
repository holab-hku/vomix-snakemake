# vomix-wrapper

A lightweight wrapper for cli integration and managing Vomix functionalities.

## Installation

```bash
cd vomix-mega
pip install .
```

## Set up conda environment

```bash
vomix activate
```

## Activate conda environment

```bash
conda activate vomix
```

## Check conda environment has been activated

```bash
conda info --envs 
```

## Usage

```bash
vomix <module> <params>
```

Example: 
```bash
vomix preprocess --outdir sample/results --datadir sample/fastq --samplelist sample/sample_list.csv
```

ctrl-C to abort

## Structure

* vomix_actions.py -> vomix actions

* vomix.py -> cli

* modules.py -> module classes 

* runModules folder -> stores the last run command for each module 

* snakemake.sh -> the script that is created and ran when running vomix <module>

## Tests

```bash
pytest tests
```