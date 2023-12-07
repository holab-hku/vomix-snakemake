rule makeblastdb_derep:
  name: "viralcontigident.py make blast db [--clustering-fast]"
  input:
    "results/viralcontigident/intermediate/scores/combined.viralcontigs.fa"
  output: 
    expand("results/viralcontigident/intermediate/derep/db.{suffix}", 
        suffix=["ndb", "nin", "not", "ntf", "nhr", "njs", "nsq", "nto"])
  params:
    outdir="results/viralcontigident/intermediate/derep/", 
    dbtype='nucl', 
    tmpdir="$TMPDIR"
  log: "logs/viralcontigident_makeblastdb.log"
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
  name: "viralcontigident.py megablast [--clustering-fast]"
  input:
    fasta="results/viralcontigident/intermediate/scores/combined.viralcontigs.fa", 
    dbcheckpoints=expand("results/viralcontigident/intermediate/derep/db.{suffix}",
                suffix=["ndb", "nin", "not", "ntf", "nhr", "njs", "nsq", "nto"])
  output:
    "results/viralcontigident/intermediate/derep/blast_out.csv"
  params:
    db="results/viralcontigident/intermediate/derep/db",
    outfmt="'6 std qlen slen'",
    maxtargetseqs=10000, 
    tmpdir="$TMPDIR"
  log: "logs/viralcontigident_megablastpairwise.log"
  conda: "../envs/checkv.yml"
  threads: 64
  resources:
    mem_mb=lambda wildcards, input, attempt: (input.size_mb) * attempt * 100
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
  name : "viralcontigident.py calculate ani [--clustering-fast]"
  input:
    "results/viralcontigident/intermediate/derep/blast_out.csv"
  output: 
    "results/viralcontigident/intermediate/derep/ani.tsv"
  params:
    script_path="workflow/scripts/viralcontigident/anicalc.py", 
    tmpdir="$TMPDIR",
  log: "logs/viralcontigident_anicalc.log"
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
  name : "viralcontigident.py cluster [--clustering-fast]"
  input:
    fa="results/viralcontigident/intermediate/scores/combined.viralcontigs.fa", 
    ani="results/viralcontigident/intermediate/derep/ani.tsv"
  output:
    tsv= "results/viralcontigident/output/derep/clusters.tsv",
    reps="results/viralcontigident/output/derep/cluster_representatives.txt"
  params:
    script="workflow/scripts/viralcontigident/aniclust.py",
    minani=config["vOTUani"], 
    targetcov=config["vOTUtargetcov"],
    querycov =config["vOTUquerycov"], 
    tmpdir="$TMPDIR", 
  log: "logs/viralcontigident_aniclust.log"
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
  name: "viralcontigident.py filter dereplicated viral contigs"
  input: 
    fna="results/viralcontigident/intermediate/scores/combined.viralcontigs.fa", 
    reps="results/viralcontigident/output/derep/cluster_representatives.txt"
  output:
    "results/viralcontigident/output/derep/combined.viralcontigs.derep.fa"
  params:
    outdir="results/viralcontigident/checkv/output",
    tmpdir="$TMPDIR/"
  log: "logs/viralcontigident_filterderep.log" 
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
