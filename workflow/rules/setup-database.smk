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
  name: "setup-database.smk VirSorter2 setup database (X.X G)"
  output: os.path.join(config['virsorter2-db'], "db.tgz")
  params:
    outdir=config['virsorter2-db'],
    tmpdir=os.path.join(tmpd, "virsorter2/db")
  log: os.path.join(logdir, "virsorter2_database.log")
  benchmark: os.path.join(benchmarks, "virsorter2_db.log")
  conda: "../envs/virsorter2.yml"
  threads: 64
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

rule chocophlan_db:
  name: "setup-database.smk HUMAnN3 chocophlan database (16.4 G)"
  output: os.path.join(config['humann-db'], "chocophlan/")
  params:
    outdir=config['humann-db'],
    tmpdir=os.path.join(tmpd, "humann/db")
  log: os.path.join(logdir, "humann_chocophlan_db.log")
  benchmark: os.path.join(benchmarks, "humann_chocophlan_db.log")
  conda: "../envs/biobakery3.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: 4000
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}
    mkdir -p {params.outdir}

    cd {params.tmpdir}
    wget --no-check-certificate http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz &> {log}
    tar -xf full_chocophlan.v201901_v31.tar.gz &> {log}
    # humann_databases --download chocophlan full {params.tmpdir}/ &> {log}

    mv {params.tmpdir}/chocophlan {params.outdir}/
    """


rule uniref_db:
  name: "setup-database.smk HUMAnN3 uniref database (20.7 G)"
  output: os.path.join(config['humann-db'], "uniref/")
  params:
    outdir=config['humann-db'],
    tmpdir=os.path.join(tmpd, "humann/db")
  log: os.path.join(logdir, "humann_uniref_db.log")
  benchmark: os.path.join(benchmarks, "humann_uniref_db.log")
  conda: "../envs/biobakery3.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: 4000
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}
    mkdir -p {params.outdir} 
 
    cd {params.tmpdir} 
    wget --no-check-certificate http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_annotated_v201901b_full.tar.gz &> {log}
    tar -xf uniref90_annotated_v201901b_full.tar.gz &> {log}
    # humann_databases --download uniref uniref90_diamond {params.tmpdir}/ &> {log}

    mv {params.tmpdir}/uniref {params.outdir}/
    """


rule utilitymap_db:
  name: "setup-database.smk HUMAnN3 utility mapping database (X.X G)"
  output: os.path.join(config['humann-db'], "utility_mapping/")
  params:
    outdir=config['humann-db'],
    tmpdir=os.path.join(tmpd, "humann/db")
  log: os.path.join(logdir, "humann_utility_db.log")
  benchmark: os.path.join(benchmarks, "humann_utility_db.log")
  conda: "../envs/biobakery3.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: 4000
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}
    mkdir -p {params.outdir} 
 
    cd {params.tmpdir} 
    wget --no-check-certificate http://huttenhower.sph.harvard.edu/humann_data/full_mapping_v201901b.tar.gz &> {log}
    tar -xf full_mapping_v201901b.tar.gz &> {log}
    # humann_databases --download utility_mapping full {params.tmpdir}/ &> {log}

    mv {params.tmpdir}/utility_mapping {params.outdir}/
    """


rule eggnog_download:
  name: "setup-database.smk eggNOG-mapper v2 Database (43.0 G)"
  output:
    os.path.join(config['eggNOG-db'], "eggnog.db"),
    os.path.join(config['eggNOG-db'], "eggnog_proteins.dmnd")
  params:
    dbdir=config['eggNOG-db'],
    outdir="workflow/database/eggNOGv2",
    tmpdir=os.path.join(tmpd, "eggNOGdb")
  conda: "../envs/eggnog-mapper.yml"
  log: os.path.join(logdir, "eggNOGv2_db.log")
  benchmark: os.path.join(benchmarks, "eggNOGv2_db.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 8 * 10**3
  shell:
    """
    rm -r {params.outdir}
    mkdir -p {params.tmpdir} {params.dbdir}
    
    download_eggnog_data.py -y --data_dir {params.tmpdir} 2> {log}
    mv {params.tmpdir}/* {params.outdir}/

    rm -rf {params.tmpdir}
    """
