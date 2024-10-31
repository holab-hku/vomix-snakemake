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
    console.print(Panel.fit("File path does not end with .fa, .fasta, or .fna", title = "Error", subtitle="Input not fasta file"))
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
      console.print(Panel.fit(f"Output directory '{outdir_p}' already exists and is not empty.", title = "Warning", subtitle="Output Directory Not Empty"))
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
    relpath("host/output/host_prediction_to_genus.csv"),
    relpath("host/output/host_prediction_to_genome.csv"),
    relpath("host/output/detailed_output_by_tool.csv")
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


rule iphop:
  name: "host.py iPHoP run"
  input:
    fasta_path
  output:
    genus=relpath("host/output/host_prediction_to_genus.csv"), 
    genome=relpath("host/output/host_prediction_to_genome.csv"), 
    detailed=relpath("host/output/detailed_output_by_tool.csv")
  params:
    outdir=relpath("host/output/"),
    dbdir=configdict['iphopdbdir'], 
    cutoff=configdict['iphopcutoff'],
    parameters=configdict['iphopparams'], 
    tmpdir=tmpd
  conda: "../envs/iphop.yml"
  log: os.path.join(logdir, "iphop.log")
  benchmark: os.path.join(benchmarks, "iphop.log")
  threads: 64
  resources: 
    mem_mb=lambda wildcards, attempt: attempt * 350 * 10**3
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    iphop predict \
        --fa_file {input} \
        --out_dir {params.tmpdir} \
        --db_dir {params.dbdir} \
        --min_score {params.cutoff} \
        --step all \
        --num_threads {threads} \
        {params.parameters} &> {log}

     mv {params.tmpdir}/Host_prediction_to_genus*.csv {output.genus}
     mv {params.tmpdir}/Host_prediction_to_genome*.csv {output.genome}
     mv {params.tmpdir}/Detailed_output_by_tool.csv file {output.detailed}

     mv {params.tmpdir}/* {params.outdir}

     rm -rf {params.tmpdir}/*

    """

