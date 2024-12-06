configdict = config['viral-annotate']

logdir = relpath("annotate/viral/logs")
tmpd = relpath("annotate/viral//tmp")
benchmarks = relpath("annotate/viral/benchmarks")

os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)

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

  outdir_p = os.path.join(cwd, relpath("taxonomy/viral/output/"))
  console.print(f"[dim]Output file will be written to the '{outdir_p}' directory.\n")

  try:
    if len(os.listdir(outdir_p)) > 0:
      console.print(Panel.fit(f"Output directory '{outdir_p}' already exists and is not empty.", title = "Warning", subtitle="Output Directory Not Empty"))
  except Exception:
    pass

  sample_id = os.path.splitext(os.path.basename(fastap))[0]

else:
  fasta_path = relpath("viralcontigident/output/combined.final.vOTUs.fa")


### MASTER RULE 
rule done_log:
  name: "viral-annotate.smk Done. removing tmp files"
  localrule: True
  input:
    os.path.join(configdict['eggNOG_db_dir'], "eggnog.db"),
    os.path.join(configdict['eggNOG_db_dir'], "eggnog_proteins.dmnd"), 
    relpath("annotate/viral/output/proteins.vOTUs.faa"), 
    relpath("annotate/viral/output/PhaVIP/final_prediction/phavip_prediction.tsv"), 
    relpath("annotate/viral/output/eggNOGv2/out.emapper.annotations")
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
rule eggnog_download:
  name: "viral-annotate.smk Download eggNOG-mapper v2 Database"
  output:
    os.path.join(configdict['eggNOG_db_dir'], "eggnog.db"), 
    os.path.join(configdict['eggNOG_db_dir'], "eggnog_proteins.dmnd")
  params:
    dbdir=configdict['eggNOG_db_dir'],
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


rule prodigalgv_taxonomy:
  name: "viral-annotate.smk prodigal-gv vTOUs [parallelized]"
  input:
    fasta_path
  output:
    relpath("annotate/viral/output/proteins.vOTUs.faa")
  params:
    script="workflow/software/prodigal-gv/parallel-prodigal-gv.py",
    outdir=relpath("annotate/viral/output"),
    tmpdir=os.path.join(tmpd, "prodigal-gv")
  conda: "../envs/prodigal-gv.yml"
  log: os.path.join(logdir, "prodigal-gv.log")
  benchmark: os.path.join(benchmarks, "prodigal-gv.log")
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 72 * 10**3
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        -i {input} \
        -a {params.tmpdir}/tmp.faa \
        -t {threads} &> {log}

    mv {params.tmpdir}/tmp.faa {output}
    rm -rf {params.tmpdir}
    """


rule eggNOGmapper:
  name: "viral-annotate.smk eggNOG-mapper v2 Run"
  input:
    faa=relpath("annotate/viral/output/proteins.vOTUs.faa"),
    db=os.path.join(configdict['eggNOG_db_dir'], "eggnog.db"),
    diamond=os.path.join(configdict['eggNOG_db_dir'], "eggnog_proteins.dmnd")
  output:
    relpath("annotate/viral/output/eggNOGv2/out.emapper.annotations")
  params:
    parameters=configdict["eggNOG_params"],
    outdir=relpath("annotate/viral/output/eggNOGv2"),
    dbdir=configdict["eggNOG_db_dir"], 
    tmpdir=os.path.join(tmpd, "eggNOGv2-mapper")
  conda: "../envs/eggnog-mapper.yml"
  log: os.path.join(logdir, "eggNOGv2-mapper.log")
  benchmark: os.path.join(benchmarks, "eggNOGv2-mapper.log")
  threads: 64
  resources: 
    mem_mb=lambda wildcards, attempt: attempt * 64 * 10**3
  shell:
    """
    rm -r {params.outdir}

    emapper.py \
        -i {input.faa} \
        --output out.emapper \
        --temp_dir {params.tmpdir}/tmp \
        --output_dir {params.tmpdir} \
        --data_dir {params.dbdir} \
        --cpu {threads} \
        {params.parameters} 2> {log}
    
    mv {params.tmpdir}/* {params.outdir}
    """




rule PhaVIP:
  name: "viral-annotate.smk PhaVIP protein annotation"
  input:
    fna=fasta_path
  output:
    relpath("annotate/viral/output/PhaVIP/final_prediction/phavip_prediction.tsv")
  params:
    parameters=configdict['PhaVIPparams'],
    dbdir=configdict['PhaVIPdb'],
    outdir=relpath("annotate/viral/output/PhaVIP"),
    tmpdir=os.path.join(tmpd, "PhaVIP")
  conda: "../envs/phabox2.yml"
  log: os.path.join(logdir, "PhaVIP.log")
  benchmark: os.path.join(benchmarks, "PhaVIP.log")
  threads: 32
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
  shell:
    """
    mkdir -p {params.tmpdir} {params.outdir}

    phabox2 --task phavip \
        --contigs {input.fna} \
        --threads {threads} \
        --outpth {params.tmpdir} \
        --dbdir {params.dbdir} \
        {params.parameters} &> {log}

    mv {params.tmpdir}/* {params.outdir}/
    rm -rf {params.tmpdir}
    """

