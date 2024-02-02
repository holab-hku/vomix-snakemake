configfile: "config/assembly.yml"
logdir = relpath("assembly/logs")
tmpd = relpath("assembly/tmp")
benchmarks = relpath("assembly/benchmarks")

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)



# MASTER RULE
rule done:
  name: "assembly.py Done. removing tmp files"
  localrule: True
  input:
    expand(relpath("assembly/samples/{assembly_id}/output/final.contigs.fa"),  assembly_id = assemblies.keys()),
    expand(relpath("assembly/reports/{summary_type}.tsv"), summary_type = ["assemblystats", "assembly_size_dist"])
  output:
    os.path.join(logdir, "done.log")
  shell:
    """
    touch {output}
    """

### RULES 

rule megahit:
  name : "assembly.py MEGAHIT assembly"
  input:
    R1s=lambda wildcards: expand(relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"]),
    R2s=lambda wildcards: expand(relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"])
  output:
    fasta=relpath("assembly/samples/{assembly_id}/output/final.contigs.fa")
  params:
    parameters=config['megahitparams'],
    minlen=config["megahit_min_contig_len"],
    outdir=relpath("assembly/samples/{assembly_id}/output"),
    interdir=relpath("assembly/samples/{assembly_id}/intermediate/megahit"),
    tmpdir=os.path.join(tmpd, "megahit")
  log: os.path.join(logdir, "megahit_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "megahit_{assembly_id}.log")
  conda: "../envs/megahit.yml"
  threads: 24
  resources:
    mem_mb = lambda wildcards, attempt, threads, input: max(attempt * input.size_mb * 3, 2000)
  shell:
    """
    rm -rf {params.tmpdir}/{wildcards.assembly_id} {params.interdir} {params.outdir}/*
    mkdir -p {params.interdir} {params.outdir} {params.tmpdir}


    megahit \
        -1 $(echo "{input.R1s}" | tr ' ' ',') \
        -2 $(echo "{input.R2s}" | tr ' ' ',') \
        --min-contig-len {params.minlen} \
        -o {params.tmpdir}/{wildcards.assembly_id} \
        -t {threads} \
        {params.parameters} &> {log} 

    mv {params.tmpdir}/{wildcards.assembly_id}/final.contigs.fa {output.fasta}
    mv {params.tmpdir}/{wildcards.assembly_id}/* {params.interdir}
    
    """


rule assembly_stats:
  name: "assembly.py aggregate assembly statistics"
  input:
    expand(relpath("assembly/samples/{assembly_id}/output/final.contigs.fa"), assembly_id = assemblies.keys())
  output:
    stats=relpath("assembly/reports/assemblystats.tsv"),
    sizedist=relpath("assembly/reports/assembly_size_dist.tsv")
  params:
    script="workflow/scripts/assembly/assemblystats.py",
    outdir=relpath("assembly/reports"),
    tmpdir=os.path.join(tmpd, "report")
  log: os.path.join(logdir, "stats.log")
  conda: "../envs/utility.yml"
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}/* 
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} -i {input} --size-dist-file {params.tmpdir}/tmp1.tsv > {params.tmpdir}/tmp2.tsv 2> {log} 
    
    mv {params.tmpdir}/tmp1.tsv {output.sizedist}
    mv {params.tmpdir}/tmp2.tsv {output.stats}
    
    """

