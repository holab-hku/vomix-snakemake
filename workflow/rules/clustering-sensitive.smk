##################
# CD-HIT CLUSTER #
##################
rule cdhit_derep:
  name : "CD-HIT --clustering-sensitive [viralcontigident.py]"
  input:
    "results/viralcontigident/intermediate/scores/combined.viralcontigs.fa"
  output:
    fasta = "results/viralcontigident/output/derep/combined.viralcontigs.derep.fa",
    cluster = "results/viralcontigident/output/derep/combined.viralcontigs.derep.fa.clstr"
  params:
    cdhitpath = config['cdhitdir'],
    cdhitparams = config['cdhitparams'],
    output_dir = "results/viralcontigident/output/derep",
    tmp_dir = "$TMPDIR",
    tmp_file = "$TMPDIR/dereplicated.viral.contigs.fa"
  log: "logs/viralcontigident_cdhitderep.log"
  benchmark: "benchmarks/viralcontigident_cdhit.log"
  threads: 32
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 72 * 10**3
  shell:
    """
    mkdir -p {params.tmp_dir} {params.output_dir}
   
    {params.cdhitpath}cd-hit -i {input} -o {params.tmp_file} -T {threads} {params.cdhitparams} &> {log}

    mv {params.tmp_dir}/dereplicated.viral.contigs.fa {output.fasta}
    mv {params.tmp_dir}/dereplicated.viralcontigs.fa.clstr {output.cluster}
    """
