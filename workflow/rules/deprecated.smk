##############################
# VIRSORTER 2 CLASSIFICATION #
##############################


#rule vs2_classify:
#  input:
#    os.path.join(config['contigdir'], "{sample_id}/final.contigs.fa")
#  output:
#    "output/viralcontigident/{sample_id}/intermediate/vs2/final-viral-score.tsv"
#  params:
#    vs2params = config['vs2params'],
#    db_dir = config['vs2db'],
#    output_dir = "output/viralcontigident/{sample_id}/intermediate/vs2/",
#    tmp_dir = "$TMPDIR/{sample_id}"
#  log: "logs/viralcontigident_{sample_id}_vs2.log"
#  benchmark: "benchmarks/viralcontigident_{sample_id}_vs2.log"
#  conda: "pipeline/envs/vs2.yml"
#  threads: 8
#  resources:
#    runtime = lambda wildcards, attempt: attempt*attempt*60
#  shell:
#    """
#    rm -rf {params.output_dir}
#    mkdir -p {params.tmp_dir}
#
#    virsorter run \
#        -i {input} \
#        -w {params.tmp_dir} \
#        --db-dir {params.db_dir} \
#        -j {threads} \
#        {params.vs2params} \
#        all &> {log}
#
#    mkdir -p {params.output_dir}
#    mv {params.tmp_dir}/* {params.output_dir}
#    """

#########################
# VIRBOT CLASSIFICATION #
#########################


#rule virbot_classify:
#  input:
#    os.path.join(config['contigdir'], "{sample_id}/final.contigs.fa")
#  output:
#    "output/viralcontigident/{sample_id}/intermediate/virbot/pos_contig_score.csv"
#  params:
#    virbotparams = config['virbotparams'],
#    output_dir = "output/viralcontigident/{sample_id}/intermediate/virbot/",
#    tmp_dir = "$TMPDIR/_virbot/{sample_id}",
#    tmp_dir_parent = "$TMPDIR/_virbot"
#  log: "logs/viralcontigident_{sample_id}_virbot.log"
#  benchmark: "benchmarks/viralcontigident_{sample_id}_virbot.log"
#  conda: "pipeline/envs/virbot.yml"
#  threads: 8
#  shell:
#    """
#    rm -rf {params.tmp_dir} # remove any old directories
#    rm -rf {params.output_dir}
#    mkdir -p {params.tmp_dir_parent}l
#    mkdir -p {params.output_dir}
#
#    python pipeline/bin/VirBot/VirBot.py \
#        --input {input} \
#        --output {params.tmp_dir} \
#        --threads {threads} \
#        {params.virbotparams} &> {log}
#
#    #mv {params.tmp_dir}/tmp {params.tmp_dir}/intermediate
#    mv -f {params.tmp_dir}/* {params.output_dir}
#    rm -r {params.tmp_dir}
#    """

#rule phabox_classify:
#  input:
#    os.path.join(config['contigdir'], "{sample_id}/final.contigs.fa")
#  output:
#    "output/phaboxout/{sample_id}/phamer/out/phamer_prediction.csv"
#  params:
#    phamerparams = config['phamerparams'],
#    db_dir = config['phamerdb'],
#    params_dir = "pipeline/params/phabox/",
#    output_dir = "output/phaboxout/{sample_id}/",
#    script_dir = "pipeline/bin/PhaBOX/scripts/",
#    tmp_dir = "$TMPDIR/{sample_id}"
#  log: "logs/viralcontigident_{sample_id}_phabox.log"
#  conda: "pipeline/envs/phabox.yml"
#  threads: 32
#  shell:
#    """
#    rm -rf {params.output_dir}
#    mkdir -p {params.output_dir}
#    mkdir -p {params.tmp_dir}
#
#    python pipeline/bin/PhaBOX/main.py \
#        --contigs {input} \
#        --threads {threads} \
#        --rootpth {params.tmp_dir} \
#        --dbdir {params.db_dir} \
#        --parampth {params.params_dir} \
#        --scriptpth {params.script_dir} \
#        {params.phamerparams} &> {log}
#
#    mv -f {params.tmp_dir}/* {params.output_dir}
#    """


#rule pseudophabox:
#  input:
#    expand("output/phaboxout/{sample_id}/phamer/out/phamer_prediction.csv", sample_id = samples.keys())
#  output:
#    "test.txt"
#  shell: "touch {output}"
