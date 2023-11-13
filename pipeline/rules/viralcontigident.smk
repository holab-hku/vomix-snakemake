localrules:
  merge_outputs, 
  filter_outputs, 
  cat_contigs





###########################
# GENOMAD CLASSIFICATION  #
###########################

rule genomad_classify:
  input:
    os.path.join(config['contigdir'], "{sample_id}/final.contigs.fa")
  output:
    "output/viralcontigident/{sample_id}/intermediate/genomad/final.contigs_summary/final.contigs_virus_summary.tsv"
  params:
    genomadparams = config['genomadparams'],
    db_dir = "pipeline/database/genomad",
    output_dir = "output/viralcontigident/{sample_id}/intermediate/genomad/",
    tmp_dir = "$TMPDIR/{sample_id}"
  log: "logs/viralcontigident_{sample_id}_genomad.log"
  benchmark: "benchmarks/viralcontigident_{sample_id}_genomad.log"
  conda: "pipeline/envs/genomad.yml"
  threads: 8
  #resources:
  #runtime = lambda wildcards, attempt: attempt*attempt*60
  shell:
    """
    mkdir -p {params.output_dir}
    rm -rf {params.tmp_dir}/*
    mkdir -p {params.tmp_dir}

    genomad end-to-end \
        --cleanup \
        {input} \
        {params.tmp_dir} \
        {params.db_dir} \
        --threads {threads} \
        {params.genomadparams} &> {log}

    mv {params.tmp_dir}/* {params.output_dir}
    """



################################
# DEEPVIRFINDER CLASSIFICATION #
################################


rule dvf_classify:
  input:
    os.path.join(config['contigdir'], "{sample_id}/final.contigs.fa")
  output:
    "output/viralcontigident/{sample_id}/intermediate/dvf/final_score.txt"
  params:
    dvfparams = config['dvfparams'], 
    model_dir = "pipeline/bin/DeepVirFinder/models/",
    output_dir = "output/viralcontigident/{sample_id}/intermediate/dvf/",
    tmp_dir = "$TMPDIR/{sample_id}"
  log: "logs/viralcontigident_{sample_id}_dvf.log"
  benchmark: "benchmarks/viralcontigident_{sample_id}_dvf.log"
  conda: "pipeline/envs/dvf.yml"
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
    os.path.join(config['contigdir'], "{sample_id}/final.contigs.fa")
  output:
    "output/viralcontigident/{sample_id}/intermediate/vs2/final-viral-score.tsv"
  params:
    vs2params = config['vs2params'],
    db_dir = config['vs2db'],
    output_dir = "output/viralcontigident/{sample_id}/intermediate/vs2/",
    tmp_dir = "$TMPDIR/{sample_id}"
  log: "logs/viralcontigident_{sample_id}_vs2.log"
  benchmark: "benchmarks/viralcontigident_{sample_id}_vs2.log"
  conda: "pipeline/envs/vs2.yml"
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
    os.path.join(config['contigdir'], "{sample_id}/final.contigs.fa")
  output:
    "output/viralcontigident/{sample_id}/intermediate/virbot/pos_contig_score.csv"
  params:
    virbotparams = config['virbotparams'],
    output_dir = "output/viralcontigident/{sample_id}/intermediate/virbot/",
    tmp_dir = "$TMPDIR/_virbot/{sample_id}",
    tmp_dir_parent = "$TMPDIR/_virbot"
  log: "logs/viralcontigident_{sample_id}_virbot.log"
  benchmark: "benchmarks/viralcontigident_{sample_id}_virbot.log"
  conda: "pipeline/envs/virbot.yml"
  threads: 8
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
    os.path.join(config['contigdir'], "{sample_id}/final.contigs.fa")
  output:
    "output/viralcontigident/{sample_id}/intermediate/phamer/out/phamer_prediction.csv"
  params:
    phamerparams = config['phamerparams'],
    db_dir = config['phamerdb'],
    params_dir = "pipeline/params/phabox/",
    output_dir = "output/viralcontigident/{sample_id}/intermediate/phamer/",
    script_dir = "pipeline/bin/PhaBOX/scripts/",
    tmp_dir = "$TMPDIR/{sample_id}"
  log: "logs/viralcontigident_{sample_id}_phamer.log"
  benchmark: "benchmarks/viralcontigident_{sample_id}_phamer.log"
  conda: "pipeline/envs/phabox.yml"
  threads: 8
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
    vs2out = "output/viralcontigident/{sample_id}/intermediate/vs2/final-viral-score.tsv",
    dvfout = "output/viralcontigident/{sample_id}/intermediate/dvf/final_score.txt", 
    virbotout = "output/viralcontigident/{sample_id}/intermediate/virbot/pos_contig_score.csv",
    phamerout = "output/viralcontigident/{sample_id}/intermediate/phamer/out/phamer_prediction.csv",
    genomadout = "output/viralcontigident/{sample_id}/intermediate/genomad/final.contigs_summary/final.contigs_virus_summary.tsv"
  output:
    "output/viralcontigident/{sample_id}/output/merged_scores.csv"
  params:
    out_dir = "output/viralcontigident/{sample_id}/output/",
    tmp_dir = "$TMPDIR/{sample_id}"
  log: "logs/viralcontigident_{sample_id}_mergeoutput.log"
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
    contig_file = os.path.join(config['contigdir'], "{sample_id}/final.contigs.fa"),
    merged_scrs = "output/viralcontigident/{sample_id}/output/merged_scores.csv"
  output:
    filtered_contigs = "output/viralcontigident/{sample_id}/output/viral.contigs.fa",
    filtered_scrs = "output/viralcontigident/{sample_id}/output/merged_scores_filtered.csv",
    positive_hits = "output/viralcontigident/{sample_id}/output/viralhits_list"
  params:
    script_path = "pipeline/src/filtercontig_scores.py",
    vs2_cutoff = config['vs2cutoff'], 
    dvf_cutoff = config['dvfcutoff'], 
    dvf_pvalmax = config['dvfpval'],
    phamer_cutoff = config['phamercutoff'], 
    phamer_pred = config['phamerpred'] 
  log: "logs/viralcontigident_{sample_id}_filtercontigs.log"
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
    expand("output/viralcontigident/{sample_id}/output/viral.contigs.fa", sample_id = sample_list)
  output: 
    "output/viralcontigident/output/combined.viralcontigs.fa"
  log: "logs/viralcontigident_catcontigs.log"
  threads: 1
  shell:
    """
    cat {input} > {output}
    """


rule cdhit_derep:
  input:
    "output/viralcontigident/output/combined.viralcontigs.fa"
  output:
    fasta = "output/viralcontigident/output/combined.viralcontigs.derep.fa",
    cluster = "output/viralcontigident/output/combined.viralcontigs.derep.fa.clstr"
  params:
    cdhitparams = config['cdhitparams'],
    coveragecutoff = config['cdhitcoverage'],
    output_dir = "output/3__viralcontigident",
    tmp_dir = "$TMPDIR",
    tmp_file = "$TMPDIR/dereplicated.viral.contigs.fa"
  log: "logs/viralcontigident_cdhitderep.log"
  benchmark: "benchmarks/viralcontigident_{sample_id}_cdhit.log"
  threads: 12
  shell:
    """
    mkdir -p {params.tmp_dir}
    mkdir -p {params.output_dir}

    cd-hit -i {input} -o {params.tmp_file} -T {threads} -c {params.coveragecutoff} {params.cdhitparams}

    mv {params.tmp_dir}/dereplicated.viral.contigs.fa {output.fasta}
    mv {params.tmp_dir}/dereplicated.viralcontigs.fa.clstr {output.cluster}
    """


rule checkv:
  input:
    "output/viralcontigident/output/combined.viralcontigs.derep.fa"
  output:
    "output/viralcontigident/output/checkv/viruses.fna"
  params:
    checkvparams= config['checkvparams'],
    output_dir = "output/viralcontigident/checkv/output",
    tmp_dir = "$TMPDIR/checkv",
    db_dir = "pipeline/database/checkv"
  log: "logs/viralcontigident_cdhitderep.log"
  benchmark: "benchmarks/viralcontigident_{sample_id}_checkv.log"
  threads: 64
  conda: "pipeline/envs/checkv.yml"
  shell:
    """
    rm -rf {params.output_dir}/*
    mkdir -p {params.output_dir}
    rm -rf {params.tmp_dir}
    mkdir -p {params.tmp_dir}

    checkv end_to_end {input} {params.tmp_dir} -d {params.db_dir} -t {threads} {params.checkvparams}

    mv {params.tmp_dir}/* {params.output_dir}/
    """
    
    

