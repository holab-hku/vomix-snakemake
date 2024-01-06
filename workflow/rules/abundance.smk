import os

configfile: "config/abundance.yml"
logdir = relpath("abundance/logs")
tmpd = relpath("abundance/tmp")
os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

methodslist = config["covermmethods"].split()
methods_c = ",".join(methodslist)

rule done:
  name: "abundance.py Done. removing tmp files"
  input:
    pseudo=expand(relpath("abundance/output/vOTU_table_{methods}.tsv"), methods = methodslist)
  output:
    os.path.join(logdir, "done.log")
  params:
    tmpdir=tmpd,
    interdir=relpath("abundance/intermediate")
  log: os.path.join(logdir, "done.log")
  shell:
    """
    rm -rf {params.tmpdir}/* {params.interdir}/*
    touch {output}
    """

rule coverm_endtoend:
  name: "abundance.py CoverM calculate abundance"
  input:
    vOTUs=relpath("viralcontigident/output/combined.final.vOTUs.fa"),
    R1=relpath("preprocess/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
    R2=relpath("preprocess/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz")
  output:
    tsv=relpath("abundance/samples/{sample_id}/output/vOTU_table.tsv"),
    bam=relpath("abundance/samples/{sample_id}/output/{sample_id}.bam")
  params:
    parameters=config["covermparams"], 
    methods=config["covermmethods"],
    outdir=relpath("abundance/samples/{sample_id}/output"),
    bamdir=os.path.join(tmpd, "coverm/{sample_id}/bam"),
    tmpdir=os.path.join(tmpd, "coverm/{sample_id}")
  log: os.path.join(logdir, "coverm_{sample_id}.log")
  conda: "../envs/coverm.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}
  
    coverm contig \
        --coupled {input.R1} {input.R2} \
        --reference {input.vOTUs} \
        --methods {params.methods} \
         --bam-file-cache-directory {params.bamdir} \
        --output-file {params.tmpdir}/tmp.tsv \
        {params.parameters} 2> {log}

    mv {params.tmpdir}/tmp.tsv {output.tsv}
    mv {params.bamdir}/*.bam {output.bam}

    rm -rf {params.tmpdir}
    """

rule coverm_merge:
  name:  "abundance.py merge vOTU tables"
  input:
    tables=expand(relpath("abundance/samples/{sample_id}/output/vOTU_table.tsv"), sample_id = samples.keys())
  output:
    otu=expand(relpath("abundance/output/vOTU_table_{methods}.tsv"), methods = methodslist)
  params:
    script="workflow/scripts/abundance/OTU_table_parse.py",
    methods=methods_c,
    outdir=relpath("abundance/output"),
    prefix="vOTU_table",
    tmpdir=os.path.join(tmpd, "merge")
  log: os.path.join(logdir, "coverm_merge.log")
  conda: "../envs/taxonomy.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}
    
    paste {input.tables} > {params.tmpdir}/tmp.tsv

    python {params.script} \
        --covermout {params.tmpdir}/tmp.tsv \
        --methods {params.methods} \
        --outdir {params.tmpdir} \
        --filename {params.prefix} &> {log}

    rm -f {params.tmpdir}/tmp.tsv
    mv {params.tmpdir}/* {params.outdir}

    rm -rf {params.tmpdir}
    """
