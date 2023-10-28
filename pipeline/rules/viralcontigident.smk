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
    "logs/3__viralcontigident_dvf_{sample_id}.log"
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
    "output/intermediate/3__viralcontigident/vs2/{sample_id}"
  params:
    vs2params = config['vs2params'],
    db_dir = config['vs2db'],
    output_dir = "output/intermediate/3__viralcontigident/vs2",
    tmp_dir = "$TMPDIR/{sample_id}"
  log:
    "logs/3__viralcontigident_vs2_{sample_id}.log"
  conda:
    "pipeline/envs/vs2.yml"
  threads: 8
  #resources:
    #runtime = lambda wildcards, attempt: attempt*attempt*60
  shell:
    """
    rm -rf {output}
    mkidr -p {params.tmp_dir}

    virsorter run \
        -i {input} \
        -w {params.tmp_dir} \
        --db-dir {params.db_dir} \
        -j {threads} \
        {params.vs2params} \
        all &> {log}

    mkdir -p {output}
    mv {params.tmp_dir}/* {output}
    """

#########################
# VIRBOT CLASSIFICATION #
#########################


rule virbot_classify:
  input:
    os.path.join(config['contigdir'], "{sample_id}.fa")
  output:
    "output/intermediate/3__viralcontigident/virbot/{sample_id}"
  log:
    "logs/3__viralcontigident_virbot_{sample_id}.log"
  conda:
    "pipeline/envs/virbot.yml"
  threads: 8
  params:
    virbotparams = config['virbotparams'],
    output_dir = "output/intermediate/3__viralcontigident/virbot/", 
    tmp_dir = "$TMPDIR/virbot/{sample_id}"
  shell:
    """
    rm -rf {params.tmp_dir} # remove any old directories
    mkdir -p {output}
    mkdir -p {params.tmp_dir}

    python pipeline/bin/VirBot/VirBot.py \
        --input {input} \
        --output {params.tmp_dir} \
        --threads {threads} \
        {params.virbotparams} &> {log}

    mv {params.tmp_dir}/tmp {params.tmp_dir}/intermediate
    mv -f {params.tmp_dir}/* {output}
    rm -r {params.tmp_dir}
    """

#########################
# PHAMER CLASSIFICATION #
#########################

rule phamer_classify:
  input:
    os.path.join(config['contigdir'], "{sample_id}.fa")
  output:
    "output/intermediate/3__viralcontigident/phamer/{sample_id}"
  log:
    "logs/3__viralcontigident_dvf_{sample_id}.log"
  conda:
    "pipeline/envs/phabox.yml"
  threads: 8
  params:
    phamerparams = config['phamerparams'],
    db_dir = config['phamerdb'],
    params_dir = "pipeline/params/phabox/",
    output_dir = "output/intermediate/3__viralcontigident/phamer/",
    script_dir = "pipeline/bin/PhaBOX/scripts/",
    tmp_dir = "$TMPDIR/{sample_id}"
  shell:
    """
    rm -r {output}
    mkdir -p {output}
    mkdir -p {params.tmp_dir}

    python pipeline/bin/PhaBOX/PhaMer_single.py \
        --contigs {input} \
        --threads {threads} \
        --rootpth {params.tmp_dir} \
        --dbdir {params.db_dir} \
        --parampth {params.params_dir} \
        --scriptpth {params.script_dir} \
        {params.phamerparams} &> {log}

    mv -f {params.tmp_dir}/* {output}
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
    "output/3__viralcontigident/{sample_id}/merged_scores.csv"
  params:
    out_dir = "output/3__viralcontigident/{sample_id}",
    tmp_dir = "$TMPDIR/{sample_id}"
  log:
    "logs/3__viralcontigident_mergeoutput.log"
  threads: 1
  shell:
    """
    mkdir -p {params.out_dir}
    python pipeline/src/viralcontigident_mergeout.py {input.vs2out} {input.dvfout} {input.virbotout} {input.phamerout} --output {params.tmp_dir}
    mv {params.tmp_dir} {output}
    """
    
    
rule filter_output:
  input:
    contig_file = "output/2__assembly/contigs/{sample_id}.fa",
    merged_scrs = "output/3__viralcontigident/{sample_id}/merged_scores.csv"
  output:
    filtered_contigs = "output/3__viralcontigident/{sample_id}/viral_contigs.fa",
    filtered_scrs = "output/3__viralcontigident/{sample_id}/merged_scores_filtered.csv",
    positive_hits = "output/3__viralcontigident/{sample_id}/viral_hits.txt"
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
  
