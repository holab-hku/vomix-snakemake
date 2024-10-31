configdict = config['host']

logdir = relpath("host/logs")
tmpd = relpath("host/tmp")
benchmarks = relpath("host/benchmarks")

os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)


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
    relpath("viralcontigident/output/combined.final.vOTUs.fa")
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

