import os 
import json 

from rich.console import Console
from rich.progress import Progress
from rich.layout import Layout
from rich.panel import Panel
console = Console()


configdict = config['prok-binning']
logdir=relpath("binning/prokaryotic/logs")
benchmarks=relpath("binning/prokaryotic/benchmarks")
tmpd=relpath("binning/prokaryotic/tmp")

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)


n_cores = config['cores'] 
assembler = config['assembler']

############################
# Single-Sample Processing #
############################

# the "contigfile" is not the nested one

if config['inputdir']!="":

  indir = cleanpath(config['inputdir'])
  console.print(f"\nconfig['inputdir'] not empty, using '{indir}' as input for viral contig identification.")
  console.print("File names without the .fa extension will be used as sample IDs.")
  cwd = os.getcwd()
  indir_path = os.path.join(cwd, indir)

  if not os.path.exists(indir):
    console.print(Panel.fit(f"The input file path '{indir}' does not exist.", title="Error", subtitle="Contig Directory"))
    sys.exit(1)

  fasta_files = [f for f in os.listdir(indir_path) if f.endswith('.fa')]
  if len(fasta_files) == 0:
    console.print(Panel.fit(f"There are no files ending with .fa in '{indir_path}', other fasta extensions are not accepted (for now) :(", title="Error", subtitle="No .fa Files"))
    sys.exit(1) 

  assembly_ids = [os.path.basename(fasta_file).rsplit(".", 1)[0] for fasta_file in fasta_files]
  wildcards_p = os.path.join(indir, "{sample_id}.fa") 
  outdir_p = os.path.join(cwd, relpath("viralcontigident/"))
  console.print(f"Creating output directory: '{outdir_p}'.\n")

else:
  wildcards_p = relpath(os.path.join("assembly", assembler, "samples/{sample_id}/output/final.contigs.fa"))
  assembly_ids = assemblies.keys()



###########################
# Multi-sample Processing #
###########################


### MASTER RULE 

rule done_log:
  name: "prok-binning.smk Done. removing tmp files"
  localrule: True
  input:
    expand(relpath("binning/prokaryotic/samples/{sample_id}/strobealign/{sample_id}.sorted.{ext}"), sample_id=samples.keys(), ext = ["bam", "bai"]),
    expand(relpath("binning/prokaryotic/assemblies/{assembly_id}/{software}/depthfile.txt"), assembly_id=assembly_ids, software=["MetaBAT2", "MaxBin2"]),
    expand(relpath("binning/prokaryotic/assemblies/{assembly_id}/MaxBin2/bins/maxbin2.summary"), assembly_id=assembly_ids), 
    expand(relpath("binning/prokaryotic/assemblies/{assembly_id}/MetaBAT2/bins/metabat2.unbinned.fa"), assembly_id=assembly_ids), 
    expand(relpath("binning/prokaryotic/assemblies/{assembly_id}/CONCOCT/concoct_clustering_merged.csv"), assembly_id=assembly_ids),
    expand(relpath("binning/prokaryotic/assemblies/{assembly_id}/finalbins_DASTool_summary.tsv"), assembly_id=assembly_ids),
    expand(relpath("binning/prokaryotic/output/no-drep/ids/{assembly_id}_MAGids.tsv"), assembly_id=assembly_ids),
    expand(relpath("binning/prokaryotic/output/no-drep/ids/unbinned/{assembly_id}_unbinned.fasta"), assembly_id=assembly_ids),
    relpath("binning/prokaryotic/output/drep/figures/Primary_clustering_dendrogram.pdf")
  output:
    os.path.join(logdir, "done.log")
  params:
    tmpdir=tmpd
  log: os.path.join(logdir, "done.log")
  shell:
    """
    rm -rf {params.tmpdir}/*
    touch {output}
    """


### RULES

rule strobealign:
  name: "prok-binning.smk strobealign SR mapping"
  input:
    R1=relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"),
    R2=relpath("preprocess/samples/{sample_id}/{sample_id}_R2.fastq.gz"),
    fasta=lambda wildcards: relpath(os.path.join("assembly", assembler, "samples", samples[wildcards.sample_id]["assembly"], "output/final.contigs.fa")),
  output:
    bam=relpath("binning/prokaryotic/samples/{sample_id}/strobealign/{sample_id}.sorted.bam")
  params:
    strobealignparams=configdict["strobealignparams"],
    tmpdir=os.path.join(tmpd, "strobealign/{sample_id}"),
  log: os.path.join(logdir, "strobealign_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "strobealign_{sample_id}.log")
  conda: "../envs/strobealign.yml"
  threads: min(4, n_cores)
  resources:
    mem_mb=lambda wildcards, attempt, input: 8 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}

    strobealign \
        -t {threads} \
        {input.fasta} \
        {input.R1} \
        {input.R2} \
        {params.strobealignparams} 2>{log} | samtools sort - -o {params.tmpdir}/tmp.bam 2> {log}
    
    mv {params.tmpdir}/tmp.bam {output}
    """

rule indexbam:
  name: "prok-binning.smk index sorted BAM"
  localrule: True
  input:
    relpath("binning/prokaryotic/samples/{sample_id}/strobealign/{sample_id}.sorted.bam")
  output:
    relpath("binning/prokaryotic/samples/{sample_id}/strobealign/{sample_id}.sorted.bai")
  log: os.path.join(logdir, "indexbam_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "indexbam_{sample_id}.log")
  conda: "../envs/strobealign.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
  shell:
    """
    samtools index {input} -o {output} 2> {log}
    """

rule binprep:
  name: "prok-binning.smk jgi summarize bam contig depths"
  input:
    bams=lambda wildcards: expand(relpath("binning/prokaryotic/samples/{sample_id}/strobealign/{sample_id}.sorted.bam"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"]),
    bais=lambda wildcards: expand(relpath("binning/prokaryotic/samples/{sample_id}/strobealign/{sample_id}.sorted.bai"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"]),
  output:
    relpath("binning/prokaryotic/assemblies/{assembly_id}/MetaBAT2/depthfile.txt"), 
  params:
    script="workflow/scripts/binning/metabat2maxbin.py",
    parameters=configdict["jgi_summarize_bam_contig_depths_params"],
    outdir=relpath("binning/prokaryotic/assemblies/{assembly_id}/MetaBAT2"),
    tmpdir=os.path.join(tmpd, "MetaBAT2/{assembly_id}")
  log: os.path.join(logdir, "MetaBAT2_prep_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "MetaBAT2_prep_{assembly_id}.log")
  conda: "../envs/metabat2.yml"
  threads: min(8, n_cores)
  resources:
    mem_mb=lambda wildcards, attempt, input: 20 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir} {params.outdir}

    jgi_summarize_bam_contig_depths \
        --outputDepth {params.tmpdir}/tmp.txt \
        {input.bams} \
        {params.parameters} 2>{log}
    
    mv {params.tmpdir}/tmp.txt {output}
    """


rule metabat2maxbin:
  name: "prok-binning.smk MetaBAT2 2 MaxBin2 coverage file"
  localrule: True
  input:
    relpath("binning/prokaryotic/assemblies/{assembly_id}/MetaBAT2/depthfile.txt")
  output:
    relpath("binning/prokaryotic/assemblies/{assembly_id}/MaxBin2/depthfile.txt")
  params:
    script="workflow/scripts/binning/metabat2maxbin.py",
    outdir=relpath("binning/prokaryotic/assemblies/{assembly_id}/MaxBin2"),
  log: os.path.join(logdir, "MaxBin2_prep_{assembly_id}.log")
  conda: "../envs/utility.yml"
  threads: 1
  shell:
    """
    mkdir -p {params.outdir}

    python {params.script} \
        -inputtxt {input} \
        -outputtxt {output} \
        -ending sorted 2> {log}

    """

rule metabat2:
  name: "prok-binning.smk MetaBAT2 binning"
  input:
    fasta=relpath(os.path.join("assembly", assembler, "samples/{assembly_id}/output/final.contigs.fa")),
    txt=relpath("binning/prokaryotic/assemblies/{assembly_id}/MetaBAT2/depthfile.txt")
  output:
    relpath("binning/prokaryotic/assemblies/{assembly_id}/MetaBAT2/bins/metabat2.unbinned.fa")
  params:
    parameters=configdict["MetaBAT2params"],
    outdir=relpath("binning/prokaryotic/assemblies/{assembly_id}/MetaBAT2"), 
    tmpdir=os.path.join(tmpd, "MetaBAT2/{assembly_id}/bins")
  log: os.path.join(logdir, "MetaBAT2_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "MetaBAT2_{assembly_id}.log")
  conda: "../envs/metabat2.yml"
  threads: min(8, n_cores)
  resources:
    mem_mb=lambda wildcards, attempt, input: 6 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir} {params.outdir}/bins

    metabat2 \
        -t {threads} \
        -i {input.fasta} \
        -a {input.txt} \
        -o {params.tmpdir}/metabat2 \
        --unbinned \
        {params.parameters} 2> {log}

    mv {params.tmpdir}/* {params.outdir}/bins
    """

rule maxbin2:
  name: "prok-binning.smk MaxBin2 binning"
  input:
    jgi=relpath("binning/prokaryotic/assemblies/{assembly_id}/MaxBin2/depthfile.txt"),
    fasta=relpath(os.path.join("assembly", assembler, "samples/{assembly_id}/output/final.contigs.fa"))
  output:
    relpath("binning/prokaryotic/assemblies/{assembly_id}/MaxBin2/bins/maxbin2.summary")
  params:
    parameters=configdict["MaxBin2params"],
    outdir=relpath("binning/prokaryotic/assemblies/{assembly_id}/MaxBin2"),
    tmpdir=os.path.join(tmpd, "MaxBin2/{assembly_id}/bins")
  log: os.path.join(logdir, "MaxBin2_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "MaxBin2_{assembly_id}.log")
  conda: "../envs/maxbin2.yml"
  threads: min(8, n_cores)
  resources:
    mem_mb=lambda wildcards, attempt, input: 4 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir} {params.outdir}/bins

    run_MaxBin.pl \
        -thread {threads} \
        -contig {input.fasta} \
        -abund {input.jgi} \
        -out {params.tmpdir}/maxbin2 \
        {params.parameters} 2> {log}

    mv {params.tmpdir}/* {params.outdir}/bins
    """

rule concoctprep:
  name: "prok-binning.smk CONCOCT prepare"
  input:
    index=lambda wildcards: expand(relpath("binning/prokaryotic/samples/{sample_id}/strobealign/{sample_id}.sorted.bai"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"]),
    bams=lambda wildcards: expand(relpath("binning/prokaryotic/samples/{sample_id}/strobealign/{sample_id}.sorted.bam"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"]),
    fasta=relpath(os.path.join("assembly", assembler, "samples/{assembly_id}/output/final.contigs.fa"))
  output:
    bed=relpath("binning/prokaryotic/assemblies/{assembly_id}/CONCOCT/final.contigs.10k.bed"),
    fa=relpath("binning/prokaryotic/assemblies/{assembly_id}/CONCOCT/final.contigs.10k.fa"), 
    tsv=relpath("binning/prokaryotic/assemblies/{assembly_id}/CONCOCT/depth.tsv"),
  params:
    outdir=relpath("binning/prokaryotic/assemblies/{assembly_id}/CONCOCT"), 
    tmpdir=os.path.join(tmpd, "CONCOCT/{assembly_id}")
  log: os.path.join(logdir, "CONCOCT_prep_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "CONCOCT_prep_{assembly_id}.log")
  conda: "../envs/concoct.yml"
  threads: min(8, n_cores)
  resources:
    mem_mb=lambda wildcards, attempt, input: 2 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir} {params.outdir}

    cut_up_fasta.py {input.fasta} \
        -c 10000 \
        -o 0 \
        --merge_last \
        -b {params.tmpdir}/tmp.bed > {params.tmpdir}/tmp.fa 2> {log}

    concoct_coverage_table.py \
        {params.tmpdir}/tmp.bed \
        {input.bams} > {params.tmpdir}/tmp.tsv

    mv {params.tmpdir}/tmp.bed {output.bed}
    mv {params.tmpdir}/tmp.fa {output.fa}
    mv {params.tmpdir}/tmp.tsv {output.tsv}
    """


rule concoct:
  name: "prok-binning.smk CONCOCT binning"
  input:
    fasta=relpath(os.path.join("assembly", assembler, "samples/{assembly_id}/output/final.contigs.fa")),
    bed=relpath("binning/prokaryotic/assemblies/{assembly_id}/CONCOCT/final.contigs.10k.bed"),
    fa=relpath("binning/prokaryotic/assemblies/{assembly_id}/CONCOCT/final.contigs.10k.fa"),       
    tsv=relpath("binning/prokaryotic/assemblies/{assembly_id}/CONCOCT/depth.tsv"),
  output:
    csv=relpath("binning/prokaryotic/assemblies/{assembly_id}/CONCOCT/concoct_clustering_merged.csv")
  params:
    parameters=configdict["CONCOCTparams"],
    outdir=relpath("binning/prokaryotic/assemblies/{assembly_id}/CONCOCT"),
    tmpdir=os.path.join(tmpd, "CONCOCT/{assembly_id}")
  log: os.path.join(logdir, "CONCOCT_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "CONCOCT_{assembly_id}.log")
  conda: "../envs/concoct.yml"
  threads: min(8, n_cores)
  resources:
    mem_mb=lambda wildcards, attempt, input: 4 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}/bins {params.outdir}

    concoct \
        -t {threads} \
        --composition_file {input.fa} \
        --coverage_file {input.tsv} \
        -b {params.tmpdir}/concoct \
        {params.parameters} 2> {log}

    merge_cutup_clustering.py \
        {params.tmpdir}/concoct_clustering_gt1000.csv > \
        {params.tmpdir}/tmp.csv 2> {log}
    
    extract_fasta_bins.py {input.fasta} \
        {params.tmpdir}/tmp.csv \
        --output_path {params.tmpdir}/bins
    
    mv {params.tmpdir}/bins {params.outdir}
    mv {params.tmpdir}/tmp.csv {output.csv}
    mv {params.tmpdir}/* {params.outdir}
    """


rule contigs2bin:
  name: "prok-binning.smk Generate contigs2bin tables"
  input:
    MetaBAT2=relpath("binning/prokaryotic/assemblies/{assembly_id}/MetaBAT2/bins/metabat2.unbinned.fa"),
    MaxBin2=relpath("binning/prokaryotic/assemblies/{assembly_id}/MaxBin2/bins/maxbin2.summary"),
    CONCOCT=relpath("binning/prokaryotic/assemblies/{assembly_id}/CONCOCT/concoct_clustering_merged.csv")
  output:
    MetaBAT2=relpath("binning/prokaryotic/assemblies/{assembly_id}/MetaBAT2/contigs2bin.tsv"),
    MaxBin2=relpath("binning/prokaryotic/assemblies/{assembly_id}/MaxBin2/contigs2bin.tsv"),
    CONCOCT=relpath("binning/prokaryotic/assemblies/{assembly_id}/CONCOCT/contigs2bin.tsv")
  params:
    MetaBAT2=relpath("binning/prokaryotic/assemblies/{assembly_id}/MetaBAT2/bins"),
    MaxBin2=relpath("binning/prokaryotic/assemblies/{assembly_id}/MaxBin2/bins"),
  log: os.path.join(logdir, "contigs2bin_{assembly_id}.log")
  conda: "../envs/dastool.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
  shell:
    """
    Fasta_to_Contig2Bin.sh -i {params.MetaBAT2}/ -e fa > {output.MetaBAT2}
    Fasta_to_Contig2Bin.sh -i {params.MaxBin2}/ -e fasta > {output.MaxBin2}
    tail -n+2 {input.CONCOCT} | tr ',' '\t' > {output.CONCOCT}
    """

rule dastool:
  name: "prok-binning.smk DAS Tool consensus binning"
  input:
    fasta=relpath(os.path.join("assembly", assembler, "samples/{assembly_id}/output/final.contigs.fa")),
    MetaBAT2=relpath("binning/prokaryotic/assemblies/{assembly_id}/MetaBAT2/contigs2bin.tsv"),
    MaxBin2=relpath("binning/prokaryotic/assemblies/{assembly_id}/MaxBin2/contigs2bin.tsv"),
    CONCOCT=relpath("binning/prokaryotic/assemblies/{assembly_id}/CONCOCT/contigs2bin.tsv")
  output:
    relpath("binning/prokaryotic/assemblies/{assembly_id}/finalbins_DASTool_summary.tsv")
  params:
    parameters=configdict["DASToolparams"], 
    outdir=relpath("binning/prokaryotic/assemblies/{assembly_id}"),
    tmpdir=os.path.join(tmpd, "DASTool/{assembly_id}")
  log: os.path.join(logdir, "DASTool_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "DASTool_{assembly_id}.log")
  conda: "../envs/dastool.yml"
  threads: min(8, n_cores)
  resources:
    mem_mb=lambda wildcards, attempt, input: 8 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}

    DAS_Tool \
        -i {input.MetaBAT2},{input.MaxBin2},{input.CONCOCT} \
        -l MetaBAT2,MaxBin2,CONCOCT \
        -c {input.fasta} \
        -o {params.tmpdir}/finalbins \
        --write_bin_evals \
        -t {threads} \
        --write_bins \
        --write_unbinned \
        {params.parameters} 2> log

    mv {params.tmpdir}/* {params.outdir}
    """


rule movemags:
  name: "prok-binning.smk Combine MAGs from assemblies"
  localrule: True
  input:
    relpath("binning/prokaryotic/assemblies/{assembly_id}/finalbins_DASTool_summary.tsv")
  output:
    tsv=relpath("binning/prokaryotic/output/no-drep/ids/{assembly_id}_MAGids.tsv"), 
    unbinned=relpath("binning/prokaryotic/output/no-drep/ids/unbinned/{assembly_id}_unbinned.fasta")
  params:
    indir=relpath("binning/prokaryotic/assemblies/{assembly_id}/finalbins_DASTool_bins"),
    outdir=relpath("binning/prokaryotic/output/no-drep")
  log: os.path.join(logdir, "move_mags_{assembly_id}.log")
  benchmark: os.path.join(logdir, "move_mags_{assembly_id}.log")
  threads: 1
  shell:
    """
    mkdir -p {params.outdir}/ids {params.outdir}/unbinned
    echo -e "ids\tDASTool_id" > {output}

    mv {params.indir}/unbinned* {output.unbinned}

    ls -v {params.indir} | grep -v "unbinned" | cat -n | \
        while read n f; do \
        echo -e "{wildcards.assembly_id}_$n.fasta\t$f" >> {output.tsv}; \
        cp "{params.indir}/$f" "{params.outdir}/{wildcards.assembly_id}_$n.fasta"; \
        done 2> {log}
    """

rule drep:
  name: "prok-binning.smk dRep final bins"
  input:
    expand(relpath("binning/prokaryotic/output/no-drep/ids/{assembly_id}_MAGids.tsv"), assembly_id=assembly_ids), 
    expand(relpath("binning/prokaryotic/output/no-drep/ids/unbinned/{assembly_id}_unbinned.fasta"), assembly_id=assembly_ids)
  output:
    relpath("binning/prokaryotic/output/drep/figures/Primary_clustering_dendrogram.pdf")
  params:
    parameters=configdict["drepparams"], 
    indir=relpath("binning/prokaryotic/output/no-drep"),
    outdir=relpath("binning/prokaryotic/output/drep"), 
    tmpdir=os.path.join(tmpd, "drep")
  log: os.path.join(logdir, "drep.log")
  benchmark: os.path.join(benchmarks, "drep.log")
  conda: "../envs/drep.yml"
  threads: min(64, n_cores)
  resources:
    mem_mb=lambda wildcards, attempt, input: 8 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir} {params.outdir}

    dRep dereplicate \
        {params.tmpdir} \
        -g {params.indir}/*.fasta \
        {params.parameters} 2> {log}

    mv {params.tmpdir}/* {params.outdir}
    """
