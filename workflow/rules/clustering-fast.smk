rule makeblastdb_derep:
  name : "viralcontigident.py make blast db [--clustering-fast]"
  input:
    "results/viralcontigident/intermediate/scores/combined.viralcontigs.fa"
  output: 
    expand("results/viralcontigident/intermediate/derep/db.{suffix}", 
        suffix = ["ndb", "nin", "not", "ntf", "nhr", "njs", "nsq", "nto"])
  params:
    output_dir = "results/viralcontigident/intermediate/derep/", 
    dbtype = 'nucl', 
    tmp_dir = "$TMPDIR",
    tmp_file_prefix = "$TMPDIR/db"
  log: "logs/viralcontigident_makeblastdb.log"
  conda: "../envs/checkv.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmp_dir}/* {params.output_dir}
    mkdir -p {params.tmp_dir} {params.output_dir}

    makeblastdb -in {input} -dbtype {params.dbtype} -out {params.tmp_file_prefix} 2> {log}

    mv {params.tmp_dir}/* {params.output_dir}
    rm -rf {params.tmp_dir}/*

    """

rule megablast_derep:
  name : "viralcontigident.py megablast [--clustering-fast]"
  input:
    fasta = "results/viralcontigident/intermediate/scores/combined.viralcontigs.fa", 
    dbcheckpoints = expand("results/viralcontigident/intermediate/derep/db.{suffix}",
                suffix = ["ndb", "nin", "not", "ntf", "nhr", "njs", "nsq", "nto"])
  output:
    "results/viralcontigident/intermediate/derep/blast_out.csv"
  params:
    db = "results/viralcontigident/intermediate/derep/db",
    outfmt = "'6 std qlen slen'",
    maxtargetseqs = 10000, 
    tmp_dir = "$TMPDIR", 
    tmp_file = "$TMPDIR/tmp.tsv"
  log: "logs/viralcontigident_megablastpairwise.log"
  conda: "../envs/checkv.yml"
  threads: 64
  resources:
    mem_mb = lambda wildcards, input, attempt: (input.size_mb) * attempt * 100
  shell: 
    """
    rm -rf {params.tmp_dir}/*
    mkdir -p {params.tmp_dir}

    blastn -query {input.fasta} \
        -db {params.db} \
        -outfmt {params.outfmt} \
        -max_target_seqs {params.maxtargetseqs} \
        -out {params.tmp_file} \
        -num_threads {threads} &> {log}
    
    mv {params.tmp_file} {output}
    rm -rf {params.tmp_dir}/*

    """
  
rule anicalc_derep:
  name : "viralcontigident.py calculate ani [--clustering-fast]"
  input:
    "results/viralcontigident/intermediate/derep/blast_out.csv"
  output: 
    "results/viralcontigident/intermediate/derep/ani.tsv"
  params:
    script_path = "workflow/scripts/viralcontigident/anicalc.py", 
    tmp_dir = "$TMPDIR",
    tmp_file = "$TMPDIR/tmp.tsv"
  log: "logs/viralcontigident_anicalc.log"
  conda: "../envs/checkv.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmp_dir}/* 
    mkdir -p {params.tmp_dir}

    python {params.script_path} \
        -i {input} \
        -o {params.tmp_file} &> {log}

    mv {params.tmp_file} {output}
    rm -rf {params.tmp_dir}/*
    """


rule aniclust_derep:
  name : "viralcontigident.py cluster [--clustering-fast]"
  input:
    fasta = "results/viralcontigident/intermediate/scores/combined.viralcontigs.fa", 
    ani = "results/viralcontigident/intermediate/derep/ani.tsv"
  output:
    tsv =  "results/viralcontigident/output/derep/clusters.tsv",
    reps = "results/viralcontigident/output/derep/cluster_representatives.txt"
  params:
    script_path = "workflow/scripts/viralcontigident/aniclust.py",
    minani = config["vOTUani"], 
    targetcov = config["vOTUtargetcov"],
    querycov  = config["vOTUquerycov"], 
    tmp_dir = "$TMPDIR", 
    tmp_file = "$TMPDIR/tmp.tsv"
  log: "logs/viralcontigident_aniclust.log"
  conda: "../envs/checkv.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmp_dir}/*
    mkdir -p {params.tmp_dir}

    python {params.script_path} \
        --fna {input.fasta} \
        --ani {input.ani} \
        --out {params.tmp_file} \
        --min_ani {params.minani} \
        --min_tcov {params.targetcov} \
        --min_qcov {params.querycov} &> {log}

    mv {params.tmp_file} {output.tsv}
    cut -f1 {output.tsv} > {output.reps}
    rm -rf {params.tmp_dir}/*
    """


rule filtercontigs_derep:
  name : "viralcontigident.py filter dereplicated viral contigs"
  input: 
    fasta = "results/viralcontigident/intermediate/scores/combined.viralcontigs.fa", 
    reps = "results/viralcontigident/output/derep/cluster_representatives.txt"
  output:
    "results/viralcontigident/output/derep/combined.viralcontigs.derep.fa"
  params:
    output_dir = "results/viralcontigident/checkv/output",
    tmp_dir = "$TMPDIR/",
    tmp_file = "$TMPDIR/tmp.fa"
  log: "logs/viralcontigident_filterderep.log" 
  conda: "../envs/utility.yml" 
  threads: 1
  shell:
    """
    rm -rf {params.tmp_dir}/*
    mkdir -p {params.tmp_dir}

    seqkit grep {input.fasta} -f {input.reps} > {params.tmp_file} 2> {log}
    mv {params.tmp_file} {output}

    rm -rf {params.tmp_dir}/*
    """
