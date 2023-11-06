import os
import sys

localrules:
  dvf_classify,
  vs2_classify, 
  virbot_classify, 
  phamer_classify, 
  merge_outputs, 
  filter_outputs

configfile: "config.yml"

# Set temporary dir
if not os.getenv("TMPDIR"):
  os.environ["TMPDIR"] = "tmp"
  os.makedirs(os.environ["TMPDIR"],exist_ok=True)


# Set wildcard constraints
wildcard_constraints:
  sample_id = "[A-Za-z0-9_\-\.]+"


#rule all:
#  input:
#    expand(config['contigdir'], "{sample_id}.fa", dataset=DATASETS)


################################
# DEEPVIRFINDER CLASSIFICATION #
################################


rule dvf_classify:
  input:
    os.path.join(config['contigdir'], "{sample_id}.fa")
  output:
    "output/intermediate/3__viralcontigident/dvf/{sample_id}/final_score.txt"
  params:
    dvfparams = config['dvfparams'], 
    model_dir = "pipeline/bin/DeepVirFinder/models/",
    output_dir = "output/intermediate/3__viralcontigident/dvf/{sample_id}",
    tmp_dir = "$TMPDIR/{sample_id}"
  log:
    "logs/3__viralcontigident_{sample_id}.log"
  conda:
    "pipeline/envs/dvf.yml"
  threads: 8
  #resources:
  #runtime = lambda wildcards, attempt: attempt*attempt*60
  shell:
    """
    mkdir -p {params.output_dir}
    mkdir -p {params.tmp_dir}

    python pipeline/bin/DeepVirFinder/dvf.py \
        -i {input} \
        -m {params.model_dir} \
        -o {params.tmp_dir} \
        {params.dvfparams} &> {log}

    mv {params.tmp_dir}/*.txt {output}
    """

##############################
# VIRSORTER 2 CLASSIFICATION #
##############################


rule vs2_classify:
  input:
    os.path.join(config['contigdir'], "{sample_id}.fa")
  output:
    "output/intermediate/3__viralcontigident/vs2/{sample_id}/final-viral-score.tsv"
  params:
    vs2params = config['vs2params'],
    db_dir = config['vs2db'],
    output_dir = "output/intermediate/3__viralcontigident/vs2/{sample_id}",
    tmp_dir = "$TMPDIR/{sample_id}"
  log:
    "logs/3__viralcontigident_{sample_id}.log"
  conda:
    "pipeline/envs/vs2.yml"
  threads: 8
  #resources:
    #runtime = lambda wildcards, attempt: attempt*attempt*60
  shell:
    """
    rm -rf {params.output_dir}
    mkdir -p {params.tmp_dir}

    virsorter run \
        -i {input} \
        -w {params.tmp_dir} \
        --db-dir {params.db_dir} \
        -j {threads} \
        {params.vs2params} \
        all &> {log}

    mkdir -p {params.output_dir}
    mv {params.tmp_dir}/* {params.output_dir}
    """

#########################
# VIRBOT CLASSIFICATION #
#########################


rule virbot_classify:
  input:
    os.path.join(config['contigdir'], "{sample_id}.fa")
  output:
    "output/intermediate/3__viralcontigident/virbot/{sample_id}/pos_contig_score.csv"
  log:
    "logs/3__viralcontigident_{sample_id}.log"
  conda:
    "pipeline/envs/virbot.yml"
  threads: 8
  params:
    virbotparams = config['virbotparams'],
    output_dir = "output/intermediate/3__viralcontigident/virbot/{sample_id}", 
    tmp_dir = "$TMPDIR/_virbot/{sample_id}",
    tmp_dir_parent = "$TMPDIR/_virbot"
  shell:
    """
    rm -rf {params.tmp_dir} # remove any old directories
    rm -rf {params.output_dir}
    mkdir -p {params.tmp_dir_parent}l
    mkdir -p {params.output_dir}

    python pipeline/bin/VirBot/VirBot.py \
        --input {input} \
        --output {params.tmp_dir} \
        --threads {threads} \
        {params.virbotparams} &> {log}

    #mv {params.tmp_dir}/tmp {params.tmp_dir}/intermediate
    mv -f {params.tmp_dir}/* {params.output_dir}
    rm -r {params.tmp_dir}
    """

#########################
# PHAMER CLASSIFICATION #
#########################

rule phamer_classify:
  input:
    os.path.join(config['contigdir'], "{sample_id}.fa")
  output:
    "output/intermediate/3__viralcontigident/phamer/{sample_id}/out/phamer_prediction.csv"
  log:
    "logs/3__viralcontigident_{sample_id}.log"
  conda:
    "pipeline/envs/phabox.yml"
  threads: 8
  params:
    phamerparams = config['phamerparams'],
    db_dir = config['phamerdb'],
    params_dir = "pipeline/params/phabox/",
    output_dir = "output/intermediate/3__viralcontigident/phamer/{sample_id}",
    script_dir = "pipeline/bin/PhaBOX/scripts/",
    tmp_dir = "tmp/{sample_id}"
  shell:
    """
    rm -rf {params.output_dir}
    mkdir -p {params.output_dir}
    mkdir -p {params.tmp_dir}

    python pipeline/bin/PhaBOX/PhaMer_single.py \
        --contigs {input} \
        --threads {threads} \
        --rootpth {params.tmp_dir} \
        --dbdir {params.db_dir} \
        --parampth {params.params_dir} \
        --scriptpth {params.script_dir} \
        {params.phamerparams} &> {log}

    mv -f {params.tmp_dir}/* {params.output_dir}
    """



#######################
# MERGED OUTPUT FILES #
#######################

rule merge_outputs:
  input:
    vs2out = "output/intermediate/3__viralcontigident/vs2/{sample_id}/final-viral-score.tsv",
    dvfout = "output/intermediate/3__viralcontigident/dvf/{sample_id}/final_score.txt", 
    virbotout = "output/intermediate/3__viralcontigident/virbot/{sample_id}/pos_contig_score.csv",
    phamerout = "output/intermediate/3__viralcontigident/phamer/{sample_id}/out/phamer_prediction.csv"
  output:
    "output/3__viralcontigident/samples/{sample_id}/merged_scores.csv"
  params:
    out_dir = "output/3__viralcontigident/samples/{sample_id}",
    tmp_dir = "tmp/{sample_id}"
  log:
    "logs/3__viralcontigident_mergeoutput.log"
  threads: 1
  shell:
    """
    mkdir -p {params.out_dir}
    mkdir -p {params.tmp_dir}
    python pipeline/src/viralcontigident_mergeout.py {input.vs2out} {input.dvfout} {input.virbotout} {input.phamerout} --output {params.tmp_dir}/tmp.csv
    mv {params.tmp_dir}/* {output}
    """
    

#########################
# FILTER ORIGINAL FASTA #
#########################

rule filter_output:
  input:
    contig_file = os.path.join(config['contigdir'], "{sample_id}.fa"),
    merged_scrs = "output/3__viralcontigident/samples/{sample_id}/merged_scores.csv"
  output:
    filtered_contigs = "output/3__viralcontigident/samples/{sample_id}/viral_contigs.fa",
    filtered_scrs = "output/3__viralcontigident/samples/{sample_id}/merged_scores_filtered.csv",
    positive_hits = "output/3__viralcontigident/samples/{sample_id}/viral_hits.txt"
  params:
    script_path = "pipeline/src/filtercontig_scores.py",
    vs2_cutoff = config['vs2cutoff'], 
    dvf_cutoff = config['dvfcutoff'], 
    dvf_pvalmax = config['dvfpval'],
    phamer_cutoff = config['phamercutoff'], 
    phamer_pred = config['phamerpred'] 
  log:
  "logs/3__viralcontigident_filtercontigs.log"
  threads: 1
  shell:
    """
    python {params.script_path} \
        --csv_path {input.merged_scrs} \
        --vs2_min_score {params.vs2_cutoff} \
        --dvf_min_score {params.dvf_cutoff} \
        --dvf_max_pval {params.dvf_pvalmax} \
        --phamer_pred {params.phamer_pred} \
        --phamer_min_score {params.phamer_cutoff} \
        --output_path {output.filtered_scrs} \
        --hitlist_path {output.positive_hits}
    
    seqtk subseq {input.contig_file} {output.positive_hits} > {output.filtered_contigs}

    """
 

##################
# CD-HIT CLUSTER #
##################

rule cat_contigs:
  input:
    expand("output/3__viralcontigident/samples/{sample_id}/viral_contigs.fa", sample_id = ["CRC_meta"])
  output: 
    "output/3__viralcontigident/combined_viralcontigs.fa"
  log:
    "logs/3__viralcontigident_catcontigs.log"
  threads: 1
  shell:
    """
    cat {input} > {output}
    """



rule derep_cluster:
  input:
   "output/3__viralcontigident/combined_viralcontigs.fa"
  output:
    "output/3__viralcontigident/output/combined_viralcontigs_derep.fa"
   #expand("output/3__viralcontigident/output/combined_viralcontigs_derep{coverage}.fa", coverage = config['cdhitcoverage'])
  params:
    cdhitparams = config['cdhitparams'],
    coveragecutoff = config['cdhitcoverage'],
    output_dir = "output/3__viralcontigident/output",
    tmp_dir = "$TMPDIR/",
    tmp_file = "$TMPDIR/dereplicated_viral_contigs.fa"
  log: 
    "logs/3__viralcontigident_cdhitderep.log"
  threads: 12
  shell:
    """
    mkdir -p {params.tmp_dir}
    mkdir -p {params.output_dir}

    cd-hit -i {input} -o {params.tmp_file} -T {threads} -c {params.coveragecutoff} {params.cdhitparams}

    mv {params.tmp_file} {output}

    
    """
