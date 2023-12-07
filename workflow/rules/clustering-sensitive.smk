##################
# CD-HIT CLUSTER #
##################
rule cdhit_derep:
  name: "CD-HIT --clustering-sensitive [viralcontigident.py]"
  input:
    "results/viralcontigident/intermediate/scores/combined.viralcontigs.fa"
  output:
    fa="results/viralcontigident/output/derep/combined.viralcontigs.derep.fa",
    clstr="results/viralcontigident/output/derep/combined.viralcontigs.derep.fa.clstr"
  params:
    cdhitpath=config['cdhitdir'],
    cdhitparams=config['cdhitparams'],
    outdir="results/viralcontigident/output/derep",
    tmpdir="$TMPDIR/cdhit",
  log: "logs/viralcontigident_cdhitderep.log"
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
