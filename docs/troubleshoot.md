# Troubleshooting Guide

## {octicon}`bug;0.85em` How to troubleshoot vomix-snakemake

If your run fails or behaves unexpectedly, follow these steps in order. This helps you identify the root cause and collect the information needed to fix the problem or report it clearly.

1. Run the command again and capture the full Snakemake output.
2. Open the job-specific log file in the results directory and the `.snakemake/log` directory.
3. Confirm the input paths, sample list, and configuration values in `config.yml` or the command line.
4. Reproduce the error using a dry run or smaller command when possible.
5. Search the repository documentation and GitHub issues for similar failures.

```{admonition} Quick tip
:class: note
When a job fails, the first error line is usually the most useful. Later lines often contain follow-on messages from the same failure.
```

## {octicon}`book;0.85em` Where to find logs

- `logs/` inside your output directory: module-specific stdout/stderr logs.
- `.snakemake/log/`: Snakemake rule execution logs and metadata.
- `.vomix/log/<timestamp>/config.json`: run metadata for the current analysis.
- Shell output printed by Snakemake on the terminal.

```{admonition} Note
:class: note
If you are using a cluster executor, scheduler logs may also be available from your cluster system (`qstat`, `squeue`, `qacct`, or the equivalent command).
```

## {octicon}`question;0.85em` Common issues and fixes

```{dropdown} Conda / environment creation fails
:open: false

- Ensure you activated the correct environment: `conda activate vomix`.
- If Snakemake cannot create environments, verify that `conda`, `mamba`, or your container engine is installed and available.
- For Apptainer mode, use `snakemake --sdm apptainer` and verify your Apptainer/Singularity setup.
```

```{dropdown} Missing input or invalid file paths
:open: false

- Check that the `fasta`, `datadir`, `samplelist`, `R1`, and `R2` values are correct and point to existing files.
- Use absolute or relative paths that are valid from the working directory where you run Snakemake.
- If a sample list CSV is used, validate the file formatting and required columns.
```

```{dropdown} Database download or path errors
:open: false

- Many modules require external databases. Ensure there is enough disk space and network access.
- If database downloads fail, rerun the appropriate setup command or module such as `module="setup-database"`.
- Check the `config.yml` keys for database paths such as `genomad-db`, `checkv-db`, `virsorter2-db`, and `vibrant-db`.
```

```{dropdown} Permission denied or filesystem issues
:open: false

- Avoid working inside directories where you do not have write permission.
- Make sure the output directory and any database directories are writable.
- If using a shared filesystem, confirm that file locking and locking proxies are supported.
```

```{dropdown} Resource limits, memory, or timeouts
:open: false

- Adjust Snakemake resources using `--resources mem_mb=<value>` or by tuning module-specific config.
- Use `splits=8` to reduce memory consumption at the cost of longer runtime.
- For cluster jobs, increase walltime, memory, or CPU resources in your scheduler submission command.
```

```{dropdown} Cluster / scheduler failures
:open: false

- Verify cluster options such as queue name, email, nodes, ppn, and memory settings in your `--cluster-generic-submit-cmd`.
- Check scheduler-specific logs for job submission failures or runtime termination.
- Confirm that the cluster login node can reach the compute nodes and that the required modules are loaded there.
```

## {octicon}`code;0.85em` Useful Snakemake commands

- `snakemake -n`: dry run, shows the planned workflow without executing.
- `snakemake --reason`: print the reason why each job is executed.
- `snakemake --rerun-incomplete`: rerun rules with incomplete outputs.
- `snakemake --printshellcmds`: show the commands executed by each rule.
- `snakemake --config module="..." ... --use-conda -j 4 --latency-wait 20`: typical command pattern for local runs.

```{admonition} If you are uncertain
:class: tip
Start with a dry run and verify your config values and file paths before executing an end-to-end workflow.
```

## {octicon}`comment-discussion;0.85em` When to raise an issue

If you cannot resolve the problem after checking logs and documentation, open an issue with the following information:

- exact command you ran
- the module name and config values you used
- the first error message from Snakemake or the failed log file
- the contents of the relevant log file(s)
- the output directory structure and `config.yml` if it is not sensitive
- the version of Snakemake, Conda, and your operating system

```{admonition} Report issue checklist
:class: note
Please include the log file and the full error message. Partial error excerpts often make it harder to identify the actual failure.
```

## {octicon}`heart;0.85em` Best practices

- Keep your workflow directory separate from analysis output.
- Use `vomix_update.sh` only to update the workflow files, not your analysis outputs.
- Re-run failed jobs with `--rerun-incomplete` before you manually delete outputs.
- Keep your Conda environment updated and install the required dependencies from `docs/install.md`.
