configfile: "config/clustering.yml"

logdir = relpath("viralcontigident/logs")
tmpd = relpath("viralcontigident/tmp")

os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

rule makeblastdb_derep:
  name: "viral-contigident.py make blast db [--clustering-fast]"
  input:
    relpath("viralcontigident/intermediate/scores/combined.viralcontigs.fa")
  output: 
    expand(relpath("viralcontigident/intermediate/derep/db.{suffix}"), 
        suffix=["ntf", "nhr", "nto"])
  params:
    outdir=relpath("viralcontigident/intermediate/derep/"), 
    dbtype='nucl', 
    tmpdir=tmpd
  log: os.path.join(logdir, "clustering/makeblastdb.log")
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
  name: "viral-contigident.py megablast [--clustering-fast]"
  input:
    fasta=relpath("viralcontigident/intermediate/scores/combined.viralcontigs.fa"), 
    dbcheckpoints=expand(relpath("viralcontigident/intermediate/derep/db.{suffix}"),
        suffix=["ntf", "nhr", "nto"])
  output:
    relpath("viralcontigident/intermediate/derep/blast_out.csv")
  params:
    db=relpath("viralcontigident/intermediate/derep/db"),
    outfmt="'6 std qlen slen'",
    maxtargetseqs=10000, 
    tmpdir=tmpd
  log: os.path.join(logdir, "clustering/megablastpairwise.log")
  conda: "../envs/checkv.yml"
  threads: 64
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
        -out {params.tmpdir}/tmp.tsv \
        -num_threads {threads} &> {log}
    
    mv {params.tmpdir}/tmp.tsv {output}
    rm -rf {params.tmpdir}/*

    """
  
rule anicalc_derep:
  name : "viral-contigident.py calculate ani [--clustering-fast]"
  input:
    relpath("viralcontigident/intermediate/derep/blast_out.csv")
  output: 
    relpath("viralcontigident/intermediate/derep/ani.tsv")
  params:
    script_path="workflow/scripts/viralcontigident/anicalc.py", 
    tmpdir=tmpd
  log: os.path.join(logdir, "clustering/anicalc.log")
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
  name : "viral-contigident.py cluster [--clustering-fast]"
  input:
    fa=relpath("viralcontigident/intermediate/scores/combined.viralcontigs.fa"), 
    ani=relpath("viralcontigident/intermediate/derep/ani.tsv")
  output:
    tsv= relpath("viralcontigident/output/derep/clusters.tsv"),
    reps=relpath("viralcontigident/output/derep/cluster_representatives.txt")
  params:
    script="workflow/scripts/viralcontigident/aniclust.py",
    minani=config["vOTUani"],
    targetcov=config["vOTUtargetcov"],
    querycov =config["vOTUquerycov"], 
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
  name: "viral-contigident.py filter dereplicated viral contigs"
  input: 
    fna=relpath("viralcontigident/intermediate/scores/combined.viralcontigs.fa"), 
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
