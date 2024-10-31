configfile: "config/community.yml"
logdir = relpath("community/bracken/logs")
tmpd = relpath("community/bracken/tmp")
benchmarks = relpath("community/bracken/benchmarks")

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True) 


rule done_log:
  name: "community.py Bracken done. Removing tmp files"
  localrule: True
  input:
  output:
  shell: 
    """
    """



rule bracken:
  name: "community.py Kraken2 "
  input: 
    R1=relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"), 
    R2=relpath("preprocess/samples/{sample_id}/{sample_id}_R2.fastq.gz"), 
  output:
    rpkm=relpath("community/bracken/samples/{sample_id}/{sample_id}.txt"),
    sam=relpath("community/bracken/samples/{sample_id}/{sample_id}_bowtie_out.txt")
  params:
    outdir=relpath("community/bracken/samples/{sample_id}"),
    bowtiedir=relpath("community/bracken/samples/{sample_id}/bowtie"),
    db=os.path.join("workflow/database/bracken"), 
    index_v=config["mpaindex_v"], 
    tmpdir=os.path.join(tmpd, "bracken/{sample_id}")
  log: os.path.join(logdir, "bracken_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "bracken_{sample_id}.log")
  threads: 16
  resources:
    mem_mb=lambda wildcards, input, attempt: attempt * 20 * 10**3
  conda: "../envs/bracken.yml"
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    bracken {input.R1},{input.R2} \
        --bowtie2db {params.db} \
        --bowtie2out {params.tmpdir}/tmp.sam \
        --input_type fastq \
        --index {params.index_v} \
        --sample_id {wildcards.sample_id} \
        -o {params.tmpdir}/tmp.txt \
        --tmp_dir {params.tmpdir} \
        --offline &> {log}

    mv {params.tmpdir}/tmp.sam {output.sam}
    mv {params.tmpdir}/tmp.txt {output.rpkm}
    """


rule metphlan_merge:
  name: "community.py merging MetaPhlAn outputs"
  localrule: True
  input:
    expand(relpath("community/bracken/samples/{sample_id}/{sample_id}.txt"), sample_id = samples.keys())
  output:
    relpath("community/bracken/output/metaphlan_out.txt")
  params:
    outdir=relpath("community/bracken/output"),
    tmpdir=os.path.join(tmpd, "bracken")
  log: os.path.join(logdir, "merge_bracken.log")
  benchmark: os.path.join(benchmarks, "merge_bracken.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, input, attempt: attempt * input.size_mb * 2
  conda: "../envs/bracken.yml"
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    merge_bracken_tables.py {input} > {params.tmpdir}/tmp.txt 

    mv {params.tmpdir}/tmp.txt {output}
    """


