configdict = config['prokaryote']
logdir = relpath("binning/prokaryotic/logs")
tmpd = relpath("binning/prokaryotic/tmp")
benchmarks = relpath("binning/prokaryotic/benchmarks")

os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True) 


rule done_log:
  name: "prokaryote.py Done. removing tmp files"
  localrule: True
  input: 
    expand(relpath("binning/prokaryotic/samples/{sample_id}/{sample_id}.sorted.bam"), sample_id = samples.keys()),
    expand(relpath("binning/prokaryotic/assemblies/{assembly_id}/proteins.faa"), assembly_id = assemblies.keys())
  output: 
    os.path.join(logdir, "done.log")
  shell:
    """
    touch {output}
    """


rule minimap2:
  name: "prokaryote.py Minimap2 short-read mapping"
  input: 
    contig=relpath("assembly/samples/{sample_id}/output/final.contigs.fa"), 
    R1=relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"), 
    R2=relpath("preprocess/samples/{sample_id}/{sample_id}_R2.fastq.gz")
  output: 
    bam=relpath("binning/prokaryotic/samples/{sample_id}/{sample_id}.sorted.bam")
  params:
    outdir=relpath("binning/prokaryotic/samples/{sample_id}"),
    tmpdir=os.path.join(tmpd, "minimap2/{sample_id}")
  benchmark: os.path.join(benchmarks, "minimap2_{sample_id}.log")
  log: os.path.join(logdir, "minimap2_{sample_id}.log")
  threads: 16
  resources: 
    mem_mb=lambda wildcards, attempt, input: attempt * 16 * 10**3
  conda: "../envs/minimap2.yml"
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    minimap2 -a -x sr {input.contig} {input.R1} {input.R2} 2>{log} \
        | samtools view -bu - \
        | samtools sort -o {params.tmpdir}/tmp.bam - 2> {log}

    mv {params.tmpdir}/tmp.bam {output.bam}
    """




rule metabat: 
  name: "prokaryote.py MetaBat2 binning"



rule prodigal_gv:
  name: "viral-binning.py prodigal-gv viral contigs"
  input:
    relpath("assembly/samples/{assembly_id}/output/final.contigs.fa")
  output:
    faa=relpath("binning/prokaryotic/intermediate/{assembly_id}/prodigal/proteins.fa")
  params:
    script="workflow/software/prodigal-gv/parallel-prodigal-gv.py",
    outdir=relpath("binning/prokaryotic/assemblies/{assembly_id}"),
    tmpdir=os.path.join(tmpd, "pyrodigal/{assembly_id}")
  conda: "../envs/prodigal-gv.yml"
  log: os.path.join(logdir, "prodigalgv_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "prodigalgv_{assembly_id}.log")
  threads: 16
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        -i {input} \
        -a {params.tmpdir}/tmp.faa \
        -t {threads} &> {log}

    mv {params.tmpdir}/tmp.faa {output.faa}
    """

rule semibin2:
  name: "prokaryote.py SemiBin2 binning"
  input:
    contigs=relpath("assembly/samples/{assembly_id}/output/final.contigs.fa"),
    proteins=relpath("binning/prokaryotic/intermediate/{assembly_id}/prodigal/proteins.fa"), 
    bams=lambda wildcards: expand(relpath("binning/prokaryotic/samples/{sample_id}/{sample_id}.sorted.bam"),
        sample_id=assemblies[wildcards.assembly_id]["sample_id"])
  output:
    relpath("binning/prokaryotic/intermediate/{assembly_id}/semibin2/recluster_bins_info.tsv")
  params:
    parameters=configdict["semibin2params"],
    outdir=relpath("binning/prokaryotic/intermediate/{assembly_id}/semibin2"),
 # log: 
 # benchmark:
  threads: 16 
  #resources:
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    SemiBin2 single_easy_bin \
        --input-fasta {input.contigs} \
        --prodigal-output-faa {input.proteins} \
        --input-bam {input.bams} \
        --tmpdir {params.tmpdir}/tmp \ 
        
    """


rule magscot:
  name: "prokaryote.py MAGScoT bin refinement"





rule metabinner:
  name: "prokaryote.py MetaBinner binning"





rule drep: 
  name: "prokaryote.py dRep dereplicate genomes"

rule checkm:
  name: "prokaryote.py CheckM quality MAGs"

rule gtdbtk: 
  name: "prokaryote.py GTDB-TK taxonomy assignment"

rule pyrogial: 
  name: "prokaryote.py Pyrodgial"

rule KEGG: 
  name: "prokaryote.py KEGG enrichmenet pathway"


