import os 

configfile: "config/clustering.yml"

logdir = relpath("viralcontigident/logs")
tmpd = relpath("viralcontigident/tmp")

os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

##################
# CD-HIT CLUSTER #
##################
rule cdhit_derep:
  name: "CD-HIT --clustering-sensitive [viralcontigident.py]"
  input:
    relpath("viralcontigident/intermediate/scores/combined.viralcontigs.fa")
  output:
    fa=relpath("/viralcontigident/output/derep/combined.viralcontigs.derep.fa"),
    clstr=relpath("/viralcontigident/output/derep/combined.viralcontigs.derep.fa.clstr)"
  params:
    cdhitpath=config['cdhitdir'],
    cdhitparams=config['cdhitparams'],
    outdir=relpath("/viralcontigident/output/derep"),
    tmpdir=os.path.join(tmpd, "cdhit")
  log: os.path.join(logdir, "clustering/cdhitderep.log")
  benchmark: "benchmarks/viralcontigident_cdhit.log"
  threads: 32
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 72 * 10**3
  shell:
    """
    mkdir -p {params.tmpdir} {params.outdir}
   
    {params.cdhitpath}/cd-hit -i {input} -o {params.tmpdir}/tmp.fa -T {threads} {params.cdhitparams} &> {log}

    mv {params.tmp_dir}/tmp.fa {output.fa}
    mv {params.tmp_dir}/tmp.fa.clstr {output.clstr}
    """
