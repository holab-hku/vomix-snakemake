import os

configfile: "config/abundance.yml"
logdir = relpath("abundance/logs")
tmpd = relpath("abundance/tmp")
os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

rule coverm_make:
  name:  "abundance.py CoverM make"
  input:
    vOTUs=relpath("viralcontigident/output/combined.final.vOTUs.fa"),
    R1=relpath("preprocess/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"), 
    R2=relpath("preprocess/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz")
  output:
    relpath("abundance/samples/{sample_id}/output/{sample_id}.bam")
  params:
    parameters=config["makeparams"],
    mapper=config["mapper"],
    outdir=relpath("abundance/samples/{sample_id}/output"),
    tmpdir=os.path.join(tmpd, "make/{sample_id}")
  log: os.path.join(logdir, "coverm_make_{sample_id}.log")
  conda: "../envs/coverm.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    coverm make \
        -1 {input.R1} \
        -2 {input.R2} \
        --reference {input.vOTUs} \
        --output-directory {params.tmpdir} \
        --mapper {params.mapper} \
        -t {threads} \
        {params.parameters} &> {log}
    
    mv {params.tmpdir}/*.bam {output}
    """

rule coverm_filter:
  name: "abundance.py CoverM filter"
  input: 
    relpath("abundance/samples/{sample_id}/output/{sample_id}.bam")
  output:
    relpath("abundance/samples/{sample_id}/output/{sample_id}.filtered.bam")
  params:
    parameters=config["filterparams"],
    outdir=relpath("abundance/samples/{sample_id}/output"), 
    tmpdir=os.path.join(tmpd, "filter/{sample_id}")
  log: os.path.join(logdir, "coverm_filter_{sample_id}.log")
  conda: "../envs/coverm.yml"
  threads: 8 
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 8 * 10**3
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}

    coverm filter \
        --bam-files {input} \
        --output-bam-files {params.tmpdir}/tmp.bam \
        --threads {threads} \
        {params.parameters} &> {log}

    mv {params.tmpdir}/tmp.bam {output}
    """

rule sort_bam:
  name: "abundance.py samtools sort"
  input:
    relpath("abundance/samples/{sample_id}/output/{sample_id}.filtered.bam")
  output:
    relpath("abundance/samples/{sample_id}/output/{sample_id}.filtered.sorted.bam")
  params:
    tmpdir=os.path.join(tmpd, "sort/{sample_id}")
  log: os.path.join(logdir, "sort_bam_{sample_id}.log")
  conda: "../envs/bowtie2.yml"
  threads: 8 
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 8 * 10**3
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}

    samtools sort --threads {threads} -o {params.tmpdir}/tmp.bam
    mv {params.tmpdir}/tmp.bam {output}
    """

rule coverm_contig:
  name: "abundance.py CoverM contig abundance"
  input:
    expand(relpath("abundance/samples/{sample_id}/output/{sample_id}.filtered.sorted.bam"), sample_id=samples.keys())
  output:
    tpm=relpath("abundance/output/vOTU_table_tmp.csv"),
    rpkm=relpath("abundance/output/vOTU_table_rpkm.csv")
  params:
    parameters=config['contigparams'],
    outdir=relpath("abundance/output/"),
    tmpdir=os.path.join(tmpd, "contig")
  log: os.path.join(logdir, "coverm.log")
  conda: "../envs/coverm.yml"
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 72 * 10**3
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    coverm contig \
        --bam-files {input} \
        --methods rpkm tpm \
        {params.parameters} > {params.tmpdir}/tmp.tsv 2> {log}

    mv {params.tmpdir}/tmp.tsv {output.rpkm}
    """
    
#rule bowtie_build:
#  name:  "abundance.py build bowtie2 index"
#  input:
#    relpath("viralcontigident/output/combined.final.vOTUs.fa")
#  output:
#    expand(relpath("abundance/intermediate/bowtie2build/final.contigs.fa.{index}.bt2l"), index = range(1,5))
#  params:
#    outdir=relpath("abundance/intermediate/bowtie2build"), 
#    prefix="final.contigs.fa",
#    tmpdir=os.path.join(tmpd, "/abundance/"
#  log: os.path.join(logdir, "bowtiebuild.log")
#  conda: "../envs/bowtie2.yml"
#  threads: 32
#  resources:
#    mem_mb=lambda wildcards, attempt: attempt * 72 * 10**3
#  shell:
#    """
#    rm -rf {params.tmpdir} {params.outdir}
#    mkdir -p {params.tmpdir} {params.outdir}
#
#    bowtie2-build \
#        --large-index \
#        --threads {threads} \
#        {input} \
#        {params.tmpdir}/{params.prefix} &> {log}
#    
#    mv {params.tmpdir}/* {params.outdir}
#    """


#rule bowtie2:
#  name: "abundance.py bowtie2 map vOTUs"
#  input:
#    R1=relpath("preprocess/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
#    R2=relpath("preprocess/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz"),
#    bowtie=expand(relpath("abundance/intermediate/bowtie2build/final.contigs.fa.{index}.bt2l"), index = range(1,5))
#  output:
#    relpath("abundance/samples/{sample_id}/output/{sample_id}.bam")
#  params:
#    parameters=config["bowtie2params"],
#    outdir=relpath("abundance/samples/{sample_id}/output/"),
#    index=relpath("abundance/intermediate/bowtie2build/final.contigs.fa"),
#    tmpdir=os.path.join(tmpd, "/abundance/{sample_id}",
#  log: os.path.join(logdir, "bowtie2_{sample_id}.log")
#  conda: "../envs/bowtie2.yml"
#  threads: 8
#  resources:
#    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
#  shell:
#    """
#    rm -rf {params.tmpdir} {params.outdir}
#    mkdir -p {params.tmpdir} {params.outdir}
#
#    bowtie2 \
#        -1 {input.R1} \
#        -2 {input.R2} \
#        -p {threads} \
#        -x {params.index} \
#        {params.parameters} 2>{log} | samtools view -h -b -o | samtools sort - -o {params.tmpdir}/tmp.bam 2> {log}
#    
#    mv {params.tmpdir}/tmp.bam {output}
#    """
