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

* vomix_actions.py -> store vomix actions (conda env set up script)

* vomix.py -> cli

* project structure: follow vomix-wrapper folder = vOMIX-MEGA folder in the original repo

## Tests

```bash
pytest tests
```