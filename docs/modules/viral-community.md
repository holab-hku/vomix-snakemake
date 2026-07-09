# `viral-community`

```{image} ../_static/figures/viral-community.svg
:width: 550
:class: no-scaled-link
:align: center
```

The `viral-community` modules estimates vOTU (viral OTU) relative abundance across your sample set. The module requires a `sample_list.csv` file, which maps the specific experimental samples to their respective paired-end read files. Utilizing the `contig` mode of `CoverM`, the module performs read mapping to quantify sequence coverage accurately.

## Quick Run

::::{tab-set}
:::{tab-item} Sample List

```bash
# Conda Run
snakemake --sdm conda --use-conda --config module="viral-community" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" -j 4 --latency-wait 20

# Apptainer (Docker Image)
snakemake --sdm apptainer --use-conda --config module="viral-community" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" -j 4 --latency-wait 20
```

:::
::::

``` {admonition} Viral Abundance Estimation
:class: Tip
Please note that, in general, viral commmunity estimations are not very robust and yet empricially tested. Each vOTU could be a viral fragment rather than representative of an entire genome, unlike prokaryotic community abundance estimation that are based on marker genes and more robust. Please proceed downstream analysis with caution. 
```

## Configuration

- **`coverm-params`:** Additional mapping or calculation flags passed to the CoverM coverage engine. (default: "--mapper minimap2-sr --min-read-percent-identity 95 --min-read-aligned-percent 75 --trim-min 10 --trim-max 90")
- **`coverm-methods`:** The calculation metric outputs selected for CoverM. Modify only for pipeline testing or debugging. (default: "tpm rpkm")

## Outputs

```bash
community/viral
├── output
    ├── vOTU_table_rpkm.tsv
    └── vOTU_table_tpm.tsv
├── samples
├── logs
├── benchmarks
└── tmp
```

The main output from this module is `vOTU_table_rpkm.tsv` and `vOTU_table_tpm.tsv`, both TSV tables summarising the estimated abundance of viral OTUS (vOTUs) identified, clustered, and annotated in the Viral Identify module. The `TPM` table is in essence a relative abundance table.

## {octicon}`book;0.85em` Troubleshooting Guide

We have specific guidelines for troubleshooting vomix-snakemake so we can help you out in your analysis journey as efficiently as possible! If you run into any unexpected errors, warnings, etc. please visit our [Troubleshooting Guide](/troubleshoot.md).

## {octicon}`bug;0.85em` Report a bug to us

Have any questions or you've found a bug during your analysis? Please don't hesitate to report it to us by making an issue on our [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/holab-hku/vomix-MEGA/issues/new).
