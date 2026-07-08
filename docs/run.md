# Run & Configuration

```{toctree}
:maxdepth: 1

modules/preprocess
modules/assembly
modules/viral-identify
modules/viral-taxonomy
modules/viral-host
modules/viral-annotate
modules/viral-community
modules/viral-benchmark
modules/viral-end-to-end
modules/prok-binning
modules/prok-annotate
modules/prok-community
modules/cluster-fast
modules/checkv-pyhmmer
modules/setup-database
```

## General Run Guide

The general workflow for running vOMIX-snakemake is to choose a module, provide input data, and pass execution options to Snakemake. Most settings are controlled through the `config/config.yml` file or by overriding them inline with `--config`.

```{admonition} Configuration File
:class: note
Use the [config.yml template](https://github.com/holab-hku/vomix-snakemake/blob/main/config/config.yml) as your reference. The YAML file must follow the expected structure, and formatting issues will trigger warnings during validation.
```

You can either:

- create a custom config file and run `snakemake --configfile config.yml`, or
- pass individual parameters directly with `--config`.

_Universal Configuration Parameters_:

- `module`: chooses the vOMIX-snakemake module to run (default: `end-to-end`)
- `workdir`: sets the working directory for the Snakefile (recommended to keep the default)
- `outdir`: selects the output directory for hierarchical results formatting (default: `./results`)
- `intermediate`: keeps large intermediate files generated during analysis (default: `False`)
- `splits`: splits data into $N$ chunks to reduce memory usage wherever possible (default: `0`)

```{admonition} Dry Run
:class: tip
Use the `-n` flag to perform a dry run first and preview the analyses that will be executed.
```

```{admonition} Snakemake Options
:class: tip
To inspect the full list of native Snakemake options, run `snakemake -h`.
```

### Command Line Format

A typical vOMIX-snakemake invocation is made up of three parts:

```bash
# 1) Snakemake command
snakemake

# 2) vOMIX-snakemake parameters
--config module="preprocess" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv"

# 3) Snakemake execution parameters
--use-conda -j 4 --latency-wait 20
```

Together they make:

```bash
snakemake --config module="preprocess" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
```

```{admonition} Input Options
:class: note
You can pass vOMIX-snakemake settings through `--config`; for a complete list, consult the [config.yml template](https://github.com/holab-hku/vomix-snakemake/blob/main/config/config.yml).
```

## General Output Structure

Use `--config outdir="results"` or set `outdir` in `config/config.yml` to define where the pipeline writes results. By default, outputs are written to `./results/` in a structured hierarchy.

```bash
{OUTDIR}
├── {MODULE NAME}
│   └── viral
│       └── output
│           ├── samples
│           ├── assemblies
│           ├── benchmarks
│           ├── logs
│           ├── tmp
│           └── reports
│   └── prokaryotic
```

_**Subdirectories**_:

1. **Output**: final cleaned and aggregated results for a module.
2. **Samples**: sample- or assembly-specific results.
3. **Assemblies**: assembly-level outputs when sample-level analyses coexist.
4. **Benchmarks**: runtime, memory, I/O, and job metadata.
5. **Logs**: standard output and error logs.
6. **Temporary**: temporary files removed after job completion.
7. **Reports**: summary reports such as assembly or preprocessing statistics.

### The `.vomix` Directory

```bash
.vomix
├── samples.json
├── assemblies.json
└── log
    ├── vomixYYYYMMDD_HHMMSS
    │   ├── sample.json
    │   ├── assemblies.json
    │   └── config.json
```

The `.vomix` directory is a hidden folder inside the results directory. It stores run-specific metadata such as the latest `sample.json`, `assemblies.json`, and a history of previous runs in the `log` directory alongside the configuration used for that run.

```{admonition} Snakemake Cache
:class: note
Snakemake also maintains its own cache and history in the `.snakemake` directory near the Snakefile. Refer to the Snakemake documentation for more details.
```

## General Input Formats

::::{tab-set}

:::{tab-item} Sample List CSV

Use a sample list for the most comprehensive workflows. This format is recommended for end-to-end analyses and maps samples to their FASTQ files or SRA accessions.

```bash
snakemake --config module="end-to-end" datadir="sample/fastq" samplelist="sample/sample_list.csv" outdir="sample/results" --use-conda -j 4 -c 4 --latency-wait 30
```

_**Options**_:

- `datadir`: points to the directory where paired-end FASTQ files are stored or where they will be downloaded.
- `samplelist`: points to a valid `sample_list.csv` file.

The `sample_list.csv` file is comma-delimited and contains the following columns:

- `sample_id`: the name of the sample
- `accession`: the SRA accession used for remote download, if needed
- `assembly`: the assembly or co-assembly name
- `R1`: the path to the forward FASTQ read
- `R2`: the path to the reverse FASTQ read

You can use `--config datadir="path/to/fastqdir"` to define where reads are stored or downloaded. The default directory is `./fastq`.

```{admonition} Sample List Examples
:class: note
The file can be used for remote downloads, co-assembly, or local FASTQ files. For local input, either provide full paths in `R1` and `R2`, or place the files in `config['datadir']` using the `<sample_id>_{1,2}.fastq.gz` naming convention.
```

:::

:::{tab-item} Single Fasta

Some modules accept a single FASTA file directly. This is useful when you want to analyze one assembly or contig set independently.

```bash
snakemake --config module="viral-identify" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results" --use-conda -j 4 --latency-wait 20
```

:::

:::{tab-item} Fasta Directory

Other modules accept a directory of FASTA files. This is convenient for batch-style analysis of multiple assemblies.

```bash
snakemake --config module="viral-identify" fastadir="sample/contigs/" outdir="sample/results" --use-conda -j 4 --latency-wait 20
```

:::

::::

```{admonition} Input Compatibility
:class: note
Not all modules support all three input types. In general, `sample_list.csv` is required for comprehensive end-to-end analyses, while `fasta` and `fastadir` allow individual modules to run independently.
```

## {octicon}`book;0.85em` Troubleshooting Guide

We have specific guidelines for troubleshooting vomix-snakemake so we can help you out in your analysis journey as efficiently as possible! If you run into any unexpected errors, warnings, etc. please visit our [Troubleshooting Guide](/troubleshoot.md).

## {octicon}`bug;0.85em` Report a bug to us

Have any questions or you've found a bug during your analysis? Please don't hesitate to report it to us by making an issue on our [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/holab-hku/vOMIX-MEGA/issues/new).
