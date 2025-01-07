logdir=relpath("annotate/prok/logs")
benchmarks=relpath("annotate/prok/benchmarks")
tmpd = relpath("annotate/prok/tmp")

email=config["email"]
nowstr=config["latest_run"]
outdir=config["outdir"] 
datadir=config["datadir"]

samples, assemblies = parse_sample_list(config["samplelist"], datadir, outdir, email, nowstr)

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)


### MASTER RULE 
rule done_log:
  name: "prok-annotate.smk Done. removing tmp files"
  localrule: True
  input:
    expand(relpath("annotate/prok/samples/{sample_id}/{sample_id}_genefamilies.tsv"), sample_id = samples.keys())
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
rule pair_fastq:
  name: "prok-annotate.smk merge paired files"
  input: 
    R1=relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"),
    R2=relpath("preprocess/samples/{sample_id}/{sample_id}_R2.fastq.gz"),
  output:
    relpath("preprocess/samples/{sample_id}/input/{sample_id}_merged.fastq")
  params:
    outdir=relpath("annotate/prok/samples/{sample_id}/input")
  conda: "../envs/seqkit-biopython.yml"
  log: os.path.join(logdir, "seqkit_concat_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "seqkit_concat_{sample_id}.log")
  threads: 1
  resources: 
    mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
  shell:
    """
    rm -r {params.outdir}
    mkdir -p {params.outdir}

    seqkit concat {input.R1} {input.R2} -f > {output}
    """


rule humann3:
  name: "prok-annotate.smk HUMAnN3 run."
  input:
    relpath("preprocess/samples/{sample_id}/input/{sample_id}_merged.fastq"),
  output:
    gene=relpath("annotate/prok/samples/{sample_id}/{sample_id}_genefamilies.tsv"),
    abund=relpath("annotate/prok/samples/{sample_id}/{sample_id}_pathabundance.tsv"),
    cov=relpath("annotate/prok/samples/{sample_id}/{sample_id}_pathcoverage.tsv")
  params:
    parameters=config["humann-params"], 
    outdir=relpath("annotate/prok/samples/{sample_id}"), 
    tmpdir=os.path.join(tmpd, "humann/{sample_id}")
  conda: "../envs/eggnog-mapper.yml"
  log: os.path.join(logdir, "HUMAnN3_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "HUMAnN3_{sample_id}.log")
  threads: 16
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 32 * 10**3
  shell:
    """
    rm -r {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}
    
    humann \
        --input {input} \
        --output {params.tmpdir} \
        --threads {threads} \
        --output-format tsv \
        {params.parameters} &> {log}

    mv {params.tmpdir}/* {params.outdir}/

    rm -rf {params.tmpdir}
    """
