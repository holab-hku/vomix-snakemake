configdict = config['viral-contigident']['clustering']
logdir = relpath("viralcontigident/logs")
tmpd = relpath("viralcontigident/tmp")
benchmarks=relpath("viralcontigident/benchmarks")


os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)


if isinstance(config['cores'], int):
 n_cores = config['cores']
else:
  console.print(Panel.fit(f"config['cores'] is not an integer: {config['cores']}, you can change the parameter in config/config.yml file", title="Error", subtitle="config['cores'] not integer"))
  sys.exit(1)


############################
# Single-Sample Processing #
############################

if config['fasta']!="":

  fastap = config['fasta']
  _, extension = os.path.splitext(fastap)

  console.print(f"\n[dim]The config['fasta'] parameter is not empty, using '{fastap}' as input.")

  if extension.lower() not in ['.fa', '.fasta', '.fna']:
    console.print(Panel.fit("File path does not end with .fa, .fasta, or .fna", title = "Error", subtitle="Input not fasta file"))
    sys.exit(1)

  cwd = os.getcwd()
  fasta_path = os.path.join(cwd, fastap)

  if not os.path.exists(fastap):
    console.print(Panel.fit("The fasta file path provided does not exist.", title="Error", subtitle="Contig File Path"))
    sys.exit(1)

  outdir_p = os.path.join(cwd, relpath("viralcontigident/output/derep/"))
  console.print(f"[dim]Output file will be written to the '{outdir_p}' directory.\n")
  
  try:
    if len(os.listdir(outdir_p)) > 0:
      console.print(Panel.fit(f"Output directory '{outdir_p}' already exists and is not empty.", title = "Warning", subtitle="Output Directory Not Empty"))
  except Exception:
    pass

  sample_id = os.path.splitext(os.path.basename(fastap))[0]

else:
  fasta_path = relpath("viralcontigident/intermediate/scores/combined.viralcontigs.fa")


### MASTER RULE

rule done_log:
  name: "clustering-fast.smk Done. removing tmp files"
  localrule: True
  input:
    expand(relpath("viralcontigident/intermediate/derep/db.{suffix}"), suffix=["ntf", "ndb"]), 
    relpath("viralcontigident/intermediate/derep/blast_out.csv"), 
    relpath("viralcontigident/intermediate/derep/ani.tsv"), 
    relpath("viralcontigident/output/derep/clusters.tsv"),
    relpath("viralcontigident/output/derep/cluster_representatives.txt"),
    relpath("viralcontigident/output/derep/combined.viralcontigs.derep.fa")
  output:
    os.path.join(logdir, "clustering-fast-done.log")
  shell: "touch {output}"



### RULES

rule makeblastdb_derep:
  name: "clustering-fast.smk make blast db [--clustering-fast]"
  input: 
    fasta_path
  output: 
    expand(relpath("viralcontigident/intermediate/derep/db.{suffix}"), suffix=["ntf", "ndb"])
  params:
    outdir=relpath("viralcontigident/intermediate/derep/"), 
    dbtype='nucl', 
    tmpdir=tmpd
  log: os.path.join(logdir, "clustering/makeblastdb.log")
  benchmark: os.path.join(benchmarks, "clustering/makeblastdb.log")
  conda: "../envs/checkv.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir}/* {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    makeblastdb -in {input} -dbtype {params.dbtype} -out {params.tmpdir}/db 2> {log}

    mv {params.tmpdir}/* {params.outdir}
    rm -rf {params.tmpdir}/*

    """

rule megablast_derep:
  name: "clustering-fast.smk megablast [--clustering-fast]"
  input:
    fasta=fasta_path, 
    dbcheckpoints=expand(relpath("viralcontigident/intermediate/derep/db.{suffix}"), suffix=["ntf", "ndb"])
  output:
    relpath("viralcontigident/intermediate/derep/blast_out.csv")
  params:
    db=relpath("viralcontigident/intermediate/derep/db"),
    outfmt="'6 std qlen slen'",
    maxtargetseqs=10000, 
    tmpdir=tmpd
  log: os.path.join(logdir, "clustering/megablastpairwise.log")
  benchmark: os.path.join(benchmarks, "clustering/megablastpairwise.log")
  conda: "../envs/checkv.yml"
  threads: min(64, n_cores)
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 72 * 10**3
  shell: 
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir}

    blastn -query {input.fasta} \
        -db {params.db} \
        -outfmt {params.outfmt} \
        -max_target_seqs {params.maxtargetseqs} \
        -out {params.tmpdir}/tmp.csv \
        -num_threads {threads} &> {log}
    
    mv {params.tmpdir}/tmp.csv {output}
    rm -rf {params.tmpdir}/*

    """
  
rule anicalc_derep:
  name : "clustering-fast.smk calculate ani [--clustering-fast]"
  input:
    relpath("viralcontigident/intermediate/derep/blast_out.csv")
  output: 
    relpath("viralcontigident/intermediate/derep/ani.tsv")
  params:
    script_path="workflow/scripts/viralcontigident/anicalc.py", 
    tmpdir=tmpd
  log: os.path.join(logdir, "clustering/anicalc.log")
  benchmark: os.path.join(benchmarks, "clustering/anicalc.log")
  conda: "../envs/checkv.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir}/* 
    mkdir -p {params.tmpdir}

    python {params.script_path} \
        -i {input} \
        -o {params.tmpdir}/tmp.tsv &> {log}

    mv {params.tmpdir}/tmp.tsv {output}
    rm -rf {params.tmpdir}/*
    """


rule aniclust_derep:
  name : "clustering-fast.smk cluster [--clustering-fast]"
  input:
    fa=fasta_path, 
    ani=relpath("viralcontigident/intermediate/derep/ani.tsv")
  output:
    tsv= relpath("viralcontigident/output/derep/clusters.tsv"),
    reps=relpath("viralcontigident/output/derep/cluster_representatives.txt")
  params:
    script="workflow/scripts/viralcontigident/aniclust.py",
    minani=configdict["vOTUani"],
    targetcov=configdict["vOTUtargetcov"],
    querycov =configdict["vOTUquerycov"], 
    tmpdir=tmpd
  log: os.path.join(logdir, "clustering/aniclust.log")
  conda: "../envs/checkv.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir}

    python {params.script} \
        --fna {input.fa} \
        --ani {input.ani} \
        --out {params.tmpdir}/tmp.tsv \
        --min_ani {params.minani} \
        --min_tcov {params.targetcov} \
        --min_qcov {params.querycov} &> {log}
    mv {params.tmpdir}/tmp.tsv {output.tsv}
    cut -f1 {output.tsv} > {output.reps}

    rm -rf {params.tmpdir}/*
    """


rule filtercontigs_derep:
  name: "clustering-fast.smk filter dereplicated viral contigs"
  input: 
    fna=fasta_path, 
    reps=relpath("viralcontigident/output/derep/cluster_representatives.txt")
  output:
    relpath("viralcontigident/output/derep/combined.viralcontigs.derep.fa")
  params:
    outdir=relpath("viralcontigident/checkv/output"),
    tmpdir=tmpd
  log: os.path.join(logdir, "clustering/filterderep.log")
  conda: "../envs/utility.yml" 
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir}

    seqkit grep {input.fna} -f {input.reps} > {params.tmpdir}/tmp.fa 2> {log}
    mv {params.tmpdir}/tmp.fa {output}

    rm -rf {params.tmpdir}/*
    """
