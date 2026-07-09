# `viral-annotate`

```{image} ../_static/figures/viral-annotate.svg
:width: 550
:class: no-scaled-link
:align: center
```

The `viral-annotate` module provides a comprehensive functional characterization of viral sequences, accepting either vOTUs derived from the `viral-identify` module or custom user-provided contig files as input. The workflow begins with protein-level gene prediction performed by `Pyrodigal-gv`, which is specifically optimized for viral genetic codes and giant virus recognition. Following prediction, the module executes a battery of annotation tools—including `eggNOG-mapper v2`, `PhaVIP`, `MetaCerberus`, and `Pharokka`—to identify protein functions, functional domains, and viral-specific features.

## Quick Run

::::{tab-set}
:::{tab-item} Sample List

```bash
# Conda Run
snakemake --sdm conda --use-conda --config module="viral-annotate" outdir="sample/results" samplelist="sample/sample_list.csv" datadir="sample/fastq" -j 4 --latency-wait 20

# Apptainer (Docker Image)
snakemake --sdm apptainer --use-conda --config module="viral-annotate" outdir="sample/results" samplelist="sample/sample_list.csv" datadir="sample/fastq" -j 4 --latency-wait 20
```

:::
:::{tab-item} Single Fasta

```bash
# Conda Run
snakemake --sdm conda --use-conda --config module="viral-annotate" outdir="sample/results" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" -j 4 --latency-wait 20

# Apptainer (Docker Image)
snakemake --sdm apptainer --use-conda --config module="viral-annotate" outdir="sample/results" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" -j 4 --latency-wait 20
```

:::
::::

``` {admonition} Dry Run
:class: Tip
Use the `-n` flag with your command to perform a dry-run first, listing all analyses that will be performed 
```

## Configuration

- **`eggNOG-params`:** Parameters for running eggNOG-mapper v2. (default: "-m diamond --hmm_evalue 0.001 --hmm_score 60 --query-cover 20 --subject-cover 20 --tax_scope auto --target_orthologs all --go_evidence non-electronic --report_orthologs")
- **`PhaVIP-params`:** Minimum contig length to filter BEFORE viral identification. (default: "")
- **`metacerberus-db`:** The directory path where the MetaCerberus database is installed or will be downloaded. (default: "database/metacerberus")
- **`metacerberus-setup-params`:** Operational configurations supplied to initialize build or index the MetaCerberus database environment. (default: "")
- **`metacerberus-params`:** Parameters for running the MetaCerberus database. (default: "--hmm ALL")
- **`pharokka-db`:** The directory path where the pharokka database is installed or will be downloaded. (default: "database/pharokka")
- **`pharokka-params`:** Additional execution parameters passed directly to the pharokka bacteriophage annotation framework. (default: "-g prodigal-gv --meta")

## Outputs

```bash
sample/results/annotate/viral
├── benchmarks
├── logs
├── output
│   ├── viral_annotate_summary.csv
│   ├── eggNOGv2
│   ├── MetaCerberus
│   ├── Pharokka
│   ├── PhaVIP
│   └── proteins.vOTUs.faa
└── tmp
```

Individual software outputs can be found in the respective subfolder. All resulting annotations are integrated into a unified, cross-referenced `.csv` file, `viral_annotate_summary.csv`, which provides a standardized summary of the functional potential of the annotated viral community.

## {octicon}`book;0.85em` Troubleshooting Guide

We have specific guidelines for troubleshooting vomix-snakemake so we can help you out in your analysis journey as efficiently as possible! If you run into any unexpected errors, warnings, etc. please visit our [Troubleshooting Guide](/troubleshoot.md).

## {octicon}`bug;0.85em` Report a bug to us

Have any questions or you've found a bug during your analysis? Please don't hesitate to report it to us by making an issue on our [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/holab-hku/vomix-MEGA/issues/new).
