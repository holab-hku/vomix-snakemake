import os

logdir="workflow/database/.logs"
benchmarks="workflow/database/.benchmarks"
tmpd="workflow/database/.tmp"

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

n_cores = config['cores']


rule done:
  name: "setup-database.smk Done."
  localrule: True
  input:
    os.path.join(config['genomad-db'], "genomad_db.source"), 
    expand("workflow/database/checkv/hmm_db/checkv_hmms/{index}.hmm", index=range(1, 81)), 
    os.path.join(config['PhaBox2-db'], "genus2hostlineage.pkl"), 
  output:
    os.path.join(benchmarks, "done.log")
  shell:
    """
    touch {output}
    """


rule genomad_db:
  name: "setup-database.smk geNomad database (1.3 G)"
  localrule: True
  output: os.path.join(config['genomad-db'], "genomad_db.source")
  params:
    outdir=config['genomad-db'],
    tmpdir=os.path.join(tmpd, "genomad/db")
  log: os.path.join(logdir, "genomad_db.log")
  benchmark: os.path.join(benchmarks, "genomad_db.log")
  conda: "../envs/genomad.yml"
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    genomad download-database {params.tmpdir} &> {log}

    mv {params.tmpdir}/genomad_db/* {params.outdir}
    """


rule checkv_db:
  name: "setup-database.smk geNomad database (7.3 G)"
  localrule: True
  output: expand("workflow/database/checkv/hmm_db/checkv_hmms/{index}.hmm", index=range(1, 81))
  params:
    outdir=config['checkv-database'],
    tmpdir=os.path.join(tmpd, "checkv/db")
  log: os.path.join(logdir, "checkv_db.log")
  benchmark: os.path.join(benchmarks, "checkv_db.log")
  conda: "../envs/checkv.yml"
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    checkv download_database {params.tmpdir} &> {log}

    mv {params.tmpdir}/checkv-db*/* {params.outdir}
    """


rule phabox2_db:
  name: "setup-database.smk PhaBox2 database (1.6 G)"
  localrule: True
  output: os.path.join(config['PhaBox2-db'], "genus2hostlineage.pkl")
  params:
    outdir=config['PhaBox2-db'],
    tmpdir=os.path.join(tmpd, "PhaBox2/db")
  log: os.path.join(logdir, "PhaBox2_db.log")
  benchmark: os.path.join(benchmarks, "PhaBox2_db.log")
  conda: "../envs/phabox2.yml"
  threads: 1
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}/phabox_db_v2
    mkdir -p {params.tmpdir} {params.outdir}

    wget -O {params.tmpdir}/phabox_db_v2.zip https://github.com/KennthShang/PhaBOX/releases/download/v2/phabox_db_v2.zip &> {log}
    unzip -o {params.tmpdir}/phabox_db_v2.zip -d {params.tmpdir} > {log}

    mv {params.tmpdir}/phabox_db_v2/* {params.outdir}
    """


rule virsorter2_db:
  name: "setup-database.smk VirSorter2 setup data"
  output: os.path.join(config['virsorter2-db'], "db.tgz")
  params:
    outdir=config['virsorter2-db'],
    tmpdir=os.path.join(tmpd, "virsorter2/db")
  log: os.path.join(logdir, "virsorter2_database.log")
  benchmark: os.path.join(benchmarks, "virsorter2_db.log")
  conda: "../envs/virsorter2.yml"
  threads: min(64, n_cores)
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: 8000
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}
    mkdir -p {params.outdir}

    virsorter setup \
        -d {params.tmpdir} \
        -j {threads} \
        -s 2> {log}

    mv {params.tmpdir}/* {params.outdir}
    """

