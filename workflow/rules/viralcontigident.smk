#localrules:
 # merge_outputs,
 #filter_outputs


###########################
# GENOMAD CLASSIFICATION  #
###########################
rule genomad_classify:
  input:
    os.path.join(config['contigdir'], "{sample_id}/output/final.contigs.fa")
  output:
    "results/viralcontigident/samples/{sample_id}/intermediate/genomad/final.contigs_summary/final.contigs_virus_summary.tsv"
  params:
    genomadparams = config['genomadparams'],
    db_dir = "workflow/database/genomad",
    output_dir = "results/viralcontigident/samples/{sample_id}/intermediate/genomad/",
    tmp_dir = "$TMPDIR/{sample_id}"
  log: "logs/viralcontigident_{sample_id}_genomad.log"
  benchmark: "benchmarks/viralcontigident_{sample_id}_genomad.log"
  conda: "../envs/genomad.yml"
  threads: 8
  resources:
    #runtime_min=lambda wildcards, attempt: attempt*attempt*120
    mem_mb=lambda wildcards, attempt: attempt * 72 * 10**3
  shell:
    """
    rm -rf {params.tmp_dir} {params.output_dir}
    mkdir -p {params.tmp_dir} {params.output_dir}

    genomad end-to-end \
        --cleanup \
        {input} \
        {params.tmp_dir} \
        {params.db_dir} \
        --threads {threads} \
        {params.genomadparams} &> {log}

    mv {params.tmp_dir}/* {params.output_dir}
    rm -rf {params.tmp_dir}
    """



################################
# DEEPVIRFINDER CLASSIFICATION #
################################


rule dvf_classify:
  input:
    os.path.join(config['contigdir'], "{sample_id}/output/final.contigs.fa")
  output:
    "results/viralcontigident/samples/{sample_id}/intermediate/dvf/final_score.txt"
  params:
    script_path = "workflow/software/DeepVirFinder/dvf.py",
    dvfparams = config['dvfparams'], 
    model_dir = "workflow/software/DeepVirFinder/models/",
    output_dir = "results/viralcontigident/samples/{sample_id}/intermediate/dvf/",
    tmp_dir = "$TMPDIR/{sample_id}"
  log: "logs/viralcontigident_{sample_id}_dvf.log"
  benchmark: "benchmarks/viralcontigident_{sample_id}_dvf.log"
  conda: "../envs/dvf.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
  shell:
    """
    mkdir -p {params.tmp_dir} {params.output_dir}

    python {params.script_path} \
        -i {input} \
        -m {params.model_dir} \
        -c {threads} \
        -o {params.tmp_dir} \
        {params.dvfparams} &> {log}

    mv {params.tmp_dir}/* {output}
    rm -rf {params.tmp_dir}
    """

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

#########################
# PHAMER CLASSIFICATION #
#########################

rule phamer_classify:
  input:
    os.path.join(config['contigdir'], "{sample_id}/output/final.contigs.fa")
  output:
    "results/viralcontigident/samples/{sample_id}/intermediate/phamer/out/phamer_prediction.csv"
  params:
    script_path = "workflow/software/PhaBOX/PhaMer_single.py",
    phamerparams = config['phamerparams'],
    db_dir = config['phamerdb'],
    params_dir = "workflow/params/phabox/",
    output_dir = "results/viralcontigident/samples/{sample_id}/intermediate/phamer/",
    script_dir = "workflow/software/PhaBOX/scripts/",
    tmp_dir = "$TMPDIR/{sample_id}"
  log: "logs/viralcontigident_{sample_id}_phamer.log"
  benchmark: "benchmarks/viralcontigident_{sample_id}_phamer.log"
  conda: "../envs/phabox.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
  shell:
    """
    rm -rf {params.output_dir}
    mkdir -p {params.tmp_dir} {params.output_dir}

    python {params.script_path} \
        --contigs {input} \
        --threads {threads} \
        --rootpth {params.tmp_dir} \
        --dbdir {params.db_dir} \
        --parampth {params.params_dir} \
        --scriptpth {params.script_dir} \
        {params.phamerparams} &> {log}

    mv -f {params.tmp_dir}/* {params.output_dir}
    rm -rf {params.tmp_dir}
    """


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


#######################
# MERGED OUTPUT FILES #
#######################

rule merge_outputs:
  input:
    genomadout = "results/viralcontigident/samples/{sample_id}/intermediate/genomad/final.contigs_summary/final.contigs_virus_summary.tsv",
    dvfout = "results/viralcontigident/samples/{sample_id}/intermediate/dvf/final_score.txt", 
    phamerout = "results/viralcontigident/samples/{sample_id}/intermediate/phamer/out/phamer_prediction.csv",
  output:
    "results/viralcontigident/samples/{sample_id}/output/merged_scores.csv"
  params:
    script_path = "workflow/scripts/viralcontigident/mergeout.py",
    genomadminlen = config['genomadminlen'],
    dvfminlen = config['dvfminlen'], 
    phamerminlen = config['dvfminlen'], 
    out_dir = "results/viralcontigident/samples/{sample_id}/output/",
    tmp_dir = "$TMPDIR/{sample_id}"
  log: "logs/viralcontigident_{sample_id}_mergeoutput.log"
  conda: "../envs/utility.yml"
  threads: 1
  shell:
    """
    mkdir -p {params.tmp_dir} {params.out_dir}

    python {params.script_path} \
        --genomadout {input.genomadout} \
        --dvfout {input.dvfout} \
        --phamerout {input.phamerout} \
        --genomadminlen {params.genomadminlen} \
        --dvfminlen {params.dvfminlen} \
        --phamerminlen {params.phamerminlen} \
        --output {params.tmp_dir}/tmp.csv &> {log}

    mv {params.tmp_dir}/* {output}
    rm -rf {params.tmp_dir}
    """
    

#########################
# FILTER ORIGINAL FASTA #
#########################

rule filter_outputs:
  input:
    contig_file = os.path.join(config['contigdir'], "{sample_id}/output/final.contigs.fa"),
    merged_scrs = "results/viralcontigident/samples/{sample_id}/output/merged_scores.csv"
  output:
    filtered_contigs = "results/viralcontigident/samples/{sample_id}/output/viral.contigs.fa",
    filtered_scrs = "results/viralcontigident/samples/{sample_id}/output/merged_scores_filtered.csv",
    positive_hits = "results/viralcontigident/samples/{sample_id}/output/viralhits_list"
  params:
    script_path = "workflow/scripts/viralcontigident/filtercontig_scores.py",
    genomad_cutoff = config['genomadcutoff'], 
    dvf_cutoff = config['dvfcutoff'], 
    dvf_pvalmax = config['dvfpval'],
    phamer_cutoff = config['phamercutoff'], 
    phamer_pred = config['phamerpred'], 
    sampleidsed = "s/^>\(.*\)$/>\\1_{sample_id}/",
    tmp_dir = "$TMPDIR/{sample_id}"
  log: "logs/viralcontigident_{sample_id}_filtercontigs.log"
  conda: "../envs/utility.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmp_dir}
    mkdir -p {params.tmp_dir}

    python {params.script_path} \
        --csv_path {input.merged_scrs} \
        --genomad_min_score {params.genomad_cutoff} \
        --dvf_min_score {params.dvf_cutoff} \
        --dvf_max_pval {params.dvf_pvalmax} \
        --phamer_pred {params.phamer_pred} \
        --phamer_min_score {params.phamer_cutoff} \
        --output_path {params.tmp_dir}/tmp.csv \
        --hitlist_path {params.tmp_dir}/tmplist

    cat {params.tmp_dir}/tmplist | uniq > {output.positive_hits}
    mv {params.tmp_dir}/tmp.csv {output.filtered_scrs}
    
    seqtk subseq {input.contig_file} {output.positive_hits} | seqtk seq -SC | sed "{params.sampleidsed}"  > {params.tmp_dir}/tmp.fa
    mv {params.tmp_dir}/tmp.fa {output.filtered_contigs}

    rm -r {params.tmp_dir}
    """
 

rule cat_contigs:
  input:
    fasta = expand("results/viralcontigident/samples/{sample_id}/output/viral.contigs.fa", sample_id = samples.keys()),
    scores = expand("results/viralcontigident/samples/{sample_id}/output/merged_scores_filtered.csv", sample_id = samples.keys())
  output: 
    fasta = "results/viralcontigident/output/combined.viralcontigs.fa",
    scores = "results/viralcontigident/output/combined.viral.scores.csv"
  params:
    script_path = "workflow/scripts/viralcontigident/mergeout_scores.py"
  log: "logs/viralcontigident_catcontigs.log"
  conda: "../envs/utility.yml"
  threads: 1
  shell:
    """
    python {params.script_path} {input.scores} > {output.scores} 2> {log}
    cat {input.fasta} > {output.fasta} 2> {log}
    """

##################
# CD-HIT CLUSTER #
##################
#rule cdhit_derep:
#  input:
#    "results/viralcontigident/output/combined.viralcontigs.fa"
#  output:
#    fasta = "results/viralcontigident/output/combined.viralcontigs.derep.fa",
#    cluster = "results/viralcontigident/output/combined.viralcontigs.derep.fa.clstr"
#  params:
#    cdhitpath = config['cdhitdir'],
#    cdhitparams = config['cdhitparams'],
#    output_dir = "results/viralcontigident/output",
#    tmp_dir = "$TMPDIR",
#    tmp_file = "$TMPDIR/dereplicated.viral.contigs.fa"
#  log: "logs/viralcontigident_cdhitderep.log"
#  benchmark: "benchmarks/viralcontigident_cdhit.log"
#  threads: 32
#  resources:
#    mem_mb = lambda wildcards, attempt: attempt * 72 * 10**3
#  shell:
#    """
#    mkdir -p {params.tmp_dir} {params.output_dir}
#    
#    {params.cdhitpath}cd-hit -i {input} -o {params.tmp_file} -T {threads} {params.cdhitparams} &> {log}
#
#    mv {params.tmp_dir}/dereplicated.viral.contigs.fa {output.fasta}
#    mv {params.tmp_dir}/dereplicated.viralcontigs.fa.clstr {output.cluster}
#    """


rule makeblastdb_derep:
  input:
    "results/viralcontigident/output/combined.viralcontigs.fa"
  output: 
    expand("results/viralcontigident/intermediate/derep/db.{suffix}", 
        suffix = ["ndb", "nin", "not", "ntf", "nhr", "njs", "nsq", "nto"])
  params:
    output_dir = "results/viralcontigident/intermediate/derep/", 
    dbtype = 'nucl', 
    tmp_dir = "$TMPDIR",
    tmp_file_prefix = "$TMPDIR/db"
  log: "logs/viralcontigident_makeblastdb.log"
  conda: "../envs/checkv.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmp_dir}/* {params.output_dir}
    mkdir -p {params.tmp_dir} {params.output_dir}

    makeblastdb -in {input} -dbtype {params.dbtype} -out {params.tmp_file_prefix} 2> {log}

    mv {params.tmp_dir}/* {params.output_dir}
    rm -rf {params.tmp_dir}/*

    """

rule megablast_derep:
  input:
    fasta = "results/viralcontigident/output/combined.viralcontigs.fa", 
    dbcheckpoints = expand("results/viralcontigident/intermediate/derep/db.{suffix}",
                suffix = ["ndb", "nin", "not", "ntf", "nhr", "njs", "nsq", "nto"])
  output:
    "results/viralcontigident/intermediate/derep/blast_out.csv"
  params:
    db = "results/viralcontigident/intermediate/derep/db",
    outfmt = "'6 std qlen slen'",
    maxtargetseqs = 10000, 
    tmp_dir = "$TMPDIR", 
    tmp_file = "$TMPDIR/tmp.tsv"
  log: "logs/viralcontigident_megablastpairwise.log"
  conda: "../envs/checkv.yml"
  threads: 64
  resources:
    mem_mb = lambda wildcards, input, attempt: (input.size_mb) * attempt * 100
  shell: 
    """
    rm -rf {params.tmp_dir}/*
    mkdir -p {params.tmp_dir}

    blastn -query {input.fasta} \
        -db {params.db} \
        -outfmt {params.outfmt} \
        -max_target_seqs {params.maxtargetseqs} \
        -out {params.tmp_file} \
        -num_threads {threads} &> {log}
    
    mv {params.tmp_file} {output}
    rm -rf {params.tmp_dir}/*

    """
  
rule anicalc_derep:
  input:
    "results/viralcontigident/intermediate/derep/blast_out.csv"
  output: 
    "results/viralcontigident/intermediate/derep/ani.tsv"
  params:
    script_path = "workflow/scripts/viralcontigident/anicalc.py", 
    tmp_dir = "$TMPDIR",
    tmp_file = "$TMPDIR/tmp.tsv"
  log: "logs/viralcontigident_anicalc.log"
  conda: "../envs/checkv.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmp_dir}/* 
    mkdir -p {params.tmp_dir}

    python {params.script_path} \
        -i {input} \
        -o {params.tmp_file} &> {log}

    mv {params.tmp_file} {output}
    rm -rf {params.tmp_dir}/*
    """


rule aniclust_derep:
  input:
    fasta = "results/viralcontigident/output/combined.viralcontigs.fa", 
    ani = "results/viralcontigident/intermediate/derep/ani.tsv"
  output:
    tsv =  "results/viralcontigident/output/derep/clusters.tsv",
    reps = "results/viralcontigident/output/derep/cluster_representatives.txt"
  params:
    script_path = "workflow/scripts/viralcontigident/aniclust.py",
    minani = config["vOTUani"], 
    targetcov = config["vOTUtargetcov"],
    querycov  = config["vOTUquerycov"], 
    tmp_dir = "$TMPDIR", 
    tmp_file = "$TMPDIR/tmp.tsv"
  log: "logs/viralcontigident_aniclust.log"
  conda: "../envs/checkv.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmp_dir}/*
    mkdir -p {params.tmp_dir}

    python {params.script_path} \
        --fna {input.fasta} \
        --ani {input.ani} \
        --out {params.tmp_file} \
        --min_ani {params.minani} \
        --min_tcov {params.targetcov} \
        --min_qcov {params.querycov} &> {log}

    mv {params.tmp_file} {output.tsv}
    cut -f1 {output.tsv} > {output.reps}
    rm -rf {params.tmp_dir}/*
    """


rule filtercontigs_derep:
  input: 
    fasta = "results/viralcontigident/output/combined.viralcontigs.fa", 
    reps = "results/viralcontigident/output/derep/cluster_representatives.txt"
  output:
    "results/viralcontigident/output/combined.viralcontigs.derep.fa"
  params:
    output_dir = "results/viralcontigident/checkv/output",
    tmp_dir = "$TMPDIR/",
    tmp_file = "$TMPDIR/tmp.fa"
  log: "logs/viralcontigident_filterderep.log" 
  conda: "../envs/utility.yml" 
  threads: 1
  shell:
    """
    rm -rf {params.tmp_dir}/*
    mkdir -p {params.tmp_dir}

    seqtk subseq {input.fasta} {input.reps} > {params.tmp_file} 2> {log}
    mv {params.tmp_file} {output}

    rm -rf {params.tmp_dir}/*
    """

rule checkv:
  input:
    "results/viralcontigident/output/combined.viralcontigs.derep.fa"
  output:
    "results/viralcontigident/output/checkv/viruses.fna"
  params:
    checkvparams= config['checkvparams'],
    output_dir = "results/viralcontigident/checkv/output",
    tmp_dir = "$TMPDIR/checkv",
    tmp_file = "$TMPDIR/tmp.fa",
    db_dir = "workflow/database/checkv"
  log: "logs/viralcontigident_cdhitderep.log"
  benchmark: "benchmarks/viralcontigident_checkv.log"
  threads: 64
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 72 * 10**3
  conda: "../envs/checkv.yml"
  shell:
    """
    rm -rf {params.tmp_dir} {params.output_dir}/*
    mkdir -p {params.tmp_dir} {params.output_dir}

    sed '/^>/ s/[_-]//g' {input} > {params.tmp_file} 2> {log}

    checkv end_to_end {params.tmp_file} {params.tmp_dir} -d {params.db_dir} -t {threads} {params.checkvparams} 2> {log}

    rm {params.tmp_file}
    mv {params.tmp_dir}/* {params.output_dir}/
    rm -rf {params.tmp_dir}
    """
    
    

