configdict = config['host']

logdir = relpath("host/logs")
tmpd = relpath("host/tmp")
benchmarks = relpath("host/benchmarks")

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
    console.print(Panel.fit("File path does not end with .fa, .fasta, or .fna", title = "Error", subtitle="Fasta Input"))
    sys.exit(1)

  cwd = os.getcwd()
  fasta_path = os.path.join(cwd, fastap)

  if not os.path.exists(fastap):
    console.print(Panel.fit("The fasta file path provided does not exist.", title="Error", subtitle="Contig File Path"))
    sys.exit(1)

  outdir_p = os.path.join(cwd, relpath("host/output/"))
  console.print(f"[dim]Output file will be written to the '{outdir_p}' directory.\n")

  try:
    if len(os.listdir(outdir_p)) > 0:
      console.print(Panel.fit(f"Output directory '{outdir_p}' already exists and is not empty.", title = "Warning", subtitle="Output Directory"))
  except Exception:
    pass

  sample_id = os.path.splitext(os.path.basename(fastap))[0]

else:
  fasta_path = relpath("viralcontigident/output/combined.final.vOTUs.fa")



### MASTER RULE 

rule done_log:
  name: "host.py Done. removing tmp files"
  localrule: True
  input:
    relpath("host/output/CHERRY/out/cherry_prediction.csv"), 
    relpath("host/output/PhaTYP/out/phatyp_prediction.csv")
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


rule CHERRY:
  name: "host.py CHERRY host prediction"
  input:
    fna=fasta_path
  output:
    relpath("host/output/CHERRY/out/cherry_prediction.csv")
  params:
    script="workflow/software/PhaBOX/Cherry_single.py",
    scriptdir="workflow/software/PhaBOX/scripts/",
    parameters=configdict['CHERRYparams'],
    paramsdir="workflow/params/phabox/",
    dbdir=configdict['CHERRYdb'],
    outdir=relpath("host/output/CHERRY"),
    tmpdir=os.path.join(tmpd, "CHERRY")
  conda: "../envs/phabox.yml"
  log: os.path.join(logdir, "CHERRY.log")
  benchmark: os.path.join(benchmarks, "CHERRY.log")
  threads: 32
  resources: 
    mem_mb=lambda wildcards, attempt: attempt * 72 * 10**3
  shell:
    """
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        --contigs {input.fna} \
        --threads {threads} \
        --rootpth {params.tmpdir} \
        --dbdir {params.dbdir} \
        --parampth {params.paramsdir} \
        --scriptpth {params.scriptdir} \
        {params.parameters} &> {log}

    mv {params.tmpdir}/* {params.outdir}/

    rm -rf {params.tmpdir}
    """

rule PhaTYP:
  name: "host.py PhaTYP lifestyle prediction"
  input:
    fna=fasta_path
  output:
    relpath("host/output/PhaTYP/out/phatyp_prediction.csv")
  params:
    script="workflow/software/PhaBOX/PhaTYP_single.py",
    scriptdir="workflow/software/PhaBOX/scripts/",
    parameters=configdict['PhaTYPparams'],
    paramsdir="workflow/params/phabox/",
    dbdir=configdict['PhaTYPdb'],
    outdir=relpath("host/output/PhaTYP"),
    tmpdir=os.path.join(tmpd, "PhaTYP")
  conda: "../envs/phabox.yml"
  log: os.path.join(logdir, "PhaTYP.log")
  benchmark: os.path.join(benchmarks, "PhaTYP.log")
  threads: 32
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
  shell:
    """
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        --contigs {input.fna} \
        --threads {threads} \
        --rootpth {params.tmpdir} \
        --dbdir {params.dbdir} \
        --parampth {params.paramsdir} \
        --scriptpth {params.scriptdir} \
        {params.parameters} &> {log}

    mv {params.tmpdir}/* {params.outdir}/

    rm -rf {params.tmpdir}
    """

