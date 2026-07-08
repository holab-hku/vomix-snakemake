# `prok-annotate`

```{image} ../_static/figures/prok-annotate.svg
:width: 550
:class: no-scaled-link
:align: center
```

Insert module summary here

## Quick Run

::::{tab-set}
:::{tab-item} Sample List

```bash
# Conda Run
snakemake --config module="setup-database" samplelist="sample/sample_list.csv" --sdm conda --use-conda -j 4 --latency-wait 20

# Apptainer (Docker Image)
snakemake --config module="setup-database" samplelist="sample/sample_list.csv" --sdm apptainer --use-conda -j 4 --latency-wait 20
```

:::
:::{tab-item} Single Fasta

```bash
# Conda Run
snakemake --config module="setup-database" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" --sdm conda --use-conda -j 4 --latency-wait 20

# Apptainer (Docker Image)
snakemake --config module="setup-database" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" --sdm apptainer --use-conda -j 4 --latency-wait 20
```

:::
:::{tab-item} Fasta Directory

```bash
# Conda Run
snakemake --config module="setup-database" fastadir="sample/contigs/" --sdm conda --use-conda -j 4 --latency-wait 20

# Apptainer (Docker Image)
snakemake --config module="setup-database" fastadir="sample/contigs/" --sdm apptainer --use-conda -j 4 --latency-wait 20
```

:::

::::

``` {admonition} Dry Run
:class: Tip
Databases are automatically setup automatically when any command that requires a database. If you would like to install all databases in one go, you may use this module. 
```

``` {admonition} Dry Run
:class: Warning
 Database sizes can be large. Run the `-n` flag to list all the databases and their estimate sizes, and make sure there is enough space on your local machine before running.
```

## Configuration

- **`hostile-index-db`:** The directory path where the Hostile database is installed or will be downloaded. (default: "database/hostile")
- **`PhaBox2-db`:** Path to PhaBox2 database for download. (default: "workflow/database/phabox_db_v2")
- **`phabox2-db-name`:** The designated database name or identifier file package required for the PhaBox2 classification tool execution. (default: "phabox_db_v2")
- **`phabox2-db-baselink`:** The primary remote server base link URL used to fetch and download resource updates for the PhaBox2 database. (default: "<https://github.com/KennthShang/PhaBOX/releases/download/v2>")
- **`genomad-db`:** Path to geNomad database for download. (default: "workflow/database/genomad")
- **`virsorter2-db`:** The directory path where the VirSorter2 database is installed or will be downloaded. (default: "database/virsorter2")
- **`vibrant-db`:** The directory path where the VIBRANT database is installed or will be downloaded. (default: "database/virsorter2")
- **`checkv-db`:** Path to CheckV database for download. (default: "workflow/database/phabox_db_v2")
- **`eggNOG-db`:** Path to eggNOG v2 database for download. (default: "workflow/database/eggNOGv2")
- **`eggNOG-db-params`:** Parameters for downloading eggNOG v2 database. (default: "")
- **`checkm2-db`:** The directory path where the CheckM2 database is installed or will be downloaded. (default: "workflow/databases/checkm2")
- **`GTDBTk-db`:** The directory path where the GTDB-Tk database is installed or will be downloaded. (default: "database/GTDB-Tk")
- **`GTDBTk-db-version`:** The reference version of the GTDB-Tk database. (default: 232)
- **`iphop-db`:** Path to iPHoP database for download. (default: "workflow/database/iphop/Aug_2023_pub_rw")
- **`iphop-db-version`:** The version identifier for the iPHoP database. (default: "iPHoP_db_Aug23_rw")
- **`iphop-db-basename`:** The primary base name of the iPHoP database. (default: "Aug_2023_pub_rw")
- **`humann-db`:** Path to HUMAnN3 databases for download. (default: "workflow/database/humann")
- **`metacerberus-db`:** The directory path where the MetaCerberus database is installed or will be downloaded. (default: "database/metacerberus")
- **`metacerberus-setup-params`:** Operational configurations supplied to initialize build or index the MetaCerberus database environment. (default: "")
- **`pharokka-db`:** The directory path where the pharokka database is installed or will be downloaded. (default: "database/pharokka")

## Outputs

```bash
{vomixdir}/workflow/database
├── checkm2
├── checkv
├── diamond
├── eggNOGv2
├── EVD
├── genomad
├── GTDB-Tk
├── hostile
├── humann
├── iphop
├── KEGG
├── kraken
├── metacerberus
├── metaphlan
├── ncbi
├── pfam
├── phabox_db_v2
├── pharokka
├── semibin
├── test
├── TIGRFAM
├── vibrant
├── viphogshmm
├── virsorter2
└── VOGDB
```

Insert output summary here

## {octicon}`book;0.85em` Troubleshooting Guide

We have specific guidelines for troubleshooting vOMIX-snakemake so we can help you out in your analysis journey as efficiently as possible! If you run into any unexpected errors, warnings, etc. please visit our [Troubleshooting Guide](/troubleshoot.md).

## {octicon}`bug;0.85em` Report a bug to us

Have any questions or you've found a bug during your analysis? Please don't hesitate to report it to us by making an issue on our [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/holab-hku/vomix-MEGA/issues/new).
