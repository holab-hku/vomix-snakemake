import os 

configdict = config['viral-identify']['clustering']
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

##################
# CD-HIT CLUSTER #
##################
rule cdhit_derep:
  name: "CD-HIT --clustering-sensitive [clustering-sensitive.smk]"
  input:
    relpath("viralcontigident/intermediate/scores/combined.viralcontigs.fa")
  output:
    fa=relpath("viralcontigident/output/derep/combined.viralcontigs.derep.fa"), 
    clstr=relpath("viralcontigident/output/derep/combined.viralcontigs.derep.fa.clstr")
  params:
    cdhitparams=configdict['cdhitparams'],
    outdir=relpath("viralcontigident/output/derep"),
    tmpdir=os.path.join(tmpd, "cdhit")
  log: os.path.join(logdir, "clustering/cdhitderep.log")
  conda: "../envs/cd-hit.yml"
  benchmark: os.path.join(benchmarks, "viralcontigident_cdhit.log")
  threads: min(32, n_cores)
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 72 * 10**3
  shell:
    """
    mkdir -p {params.tmpdir} {params.outdir}
    
    cd-hit -i {input} -o {params.tmpdir}/tmp.fa -T {threads} {params.cdhitparams} &> {log}

    mv {params.tmpdir}/tmp.fa {output.fa}
    mv {params.tmpdir}/tmp.fa.clstr {output.clstr}
    """
