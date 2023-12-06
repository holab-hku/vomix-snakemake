#localrules:
 # merge_outputs,
 #filter_outputs


###########################
# GENOMAD CLASSIFICATION  #
###########################

rule genomad_classify:
  name : "viralcontigident.py geNomad classify" 
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
    rm -rf {params.tmp_dir} {params.output_dir} 2> {log}
    mkdir -p {params.tmp_dir} {params.output_dir} 2> {log}

    genomad end-to-end \
        --cleanup \
        {input} \
        {params.tmp_dir} \
        {params.db_dir} \
        --threads {threads} \
        {params.genomadparams} &> {log}

    mv {params.tmp_dir}/* {params.output_dir} 2> {log}
    rm -rf {params.tmp_dir} 2> {log}
    """



################################
# DEEPVIRFINDER CLASSIFICATION #
################################

rule dvf_classify:
  name : "viralcontigident.py DeepVirFinder classify"
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


#########################
# PHAMER CLASSIFICATION #
#########################

rule phamer_classify:
  name : "viralcontigident.py PhaMer classify"
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



#######################
# MERGED OUTPUT FILES #
#######################

rule merge_outputs:
  name : "viralcontigident.py merge classification outputs"
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
  name : "viralcontigident.py filter viral contigs"
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
    
    seqkit grep {input.contig_file} -f {output.positive_hits} | seqkit replace -p  "\s.*" -r "" |seqkit replace -p $ -r _{wildcards.sample_id}  > {params.tmp_dir}/tmp.fa 2> {log}
    mv {params.tmp_dir}/tmp.fa {output.filtered_contigs}

    rm -r {params.tmp_dir}
    """
 

rule cat_contigs:
  name : "viralcontigident.py combine viral contigs"
  input:
    fasta = expand("results/viralcontigident/samples/{sample_id}/output/viral.contigs.fa", sample_id = samples.keys()),
    scores = expand("results/viralcontigident/samples/{sample_id}/output/merged_scores_filtered.csv", sample_id = samples.keys())
  output: 
    fasta = "results/viralcontigident/intermediate/scores/combined.viralcontigs.fa",
    scores = "results/viralcontigident/intermediate/scores/combined_viral_scores.csv"
  params:
    script_path = "workflow/scripts/viralcontigident/mergeout_scores.py", 
    tmp_dir = "$TMPDIR/viralcontigident"
  log: "logs/viralcontigident_catcontigs.log"
  conda: "../envs/utility.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmp_dir}
    mkdir -p {params.tmp_dir}

    python {params.script_path} {input.scores} > {params.tmp_dir}/tmp.csv 2> {log}
    cat {input.fasta} > {params.tmp_dir}/tmp.fa 2> {log}

    mv {params.tmp_dir}/tmp.csv {output.scores}
    mv {params.tmp_dir}/tmp.fa {output.fasta}

    rm -rf {params.tmp_dir}
    """



rule checkv:
  name : "viralcontigident.py CheckV dereplicated contigs"
  input:
    "results/viralcontigident/output/derep/combined.viralcontigs.derep.fa"
  output:
    "results/viralcontigident/output/checkv/viruses.fna",
    "results/viralcontigident/output/checkv/proviruses.fna", 
    "results/viralcontigident/output/checkv/quality_summary.tsv"
  params:
    checkvparams= config['checkvparams'],
    output_dir = "results/viralcontigident/output/checkv",
    tmp_dir = "$TMPDIR/checkv",
    db_dir = "workflow/database/checkv"
  log: "logs/viralcontigident_checkv.log"
  benchmark: "benchmarks/viralcontigident_checkv.log"
  threads: 64
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 72 * 10**3
  conda: "../envs/checkv.yml"
  shell:
    """
    rm -rf {params.tmp_dir} {params.output_dir}/*
    mkdir -p {params.tmp_dir} {params.output_dir}

    checkv end_to_end {input} {params.tmp_dir} -d {params.db_dir} -t {threads} {params.checkvparams} 2> {log}

    mv {params.tmp_dir}/* {params.output_dir}
    rm -rf {params.tmp_dir}"""


rule combine_classifications:
  name: "viralcontigident.py combine derepped classification results"
  input:
    checkv_out = "results/viralcontigident/output/checkv/quality_summary.tsv",
    classify_out = "results/viralcontigident/intermediate/scores/combined_viral_scores.csv"
  output:
    "results/viralcontigident/output/checkv/combined_classification_results.csv"
  params:
    script_path = "workflow/scripts/viralcontigident/combineclassify.py",
    tmp_dir = "$TMPDIR"
  log: "logs/viralcontigident_combineclassification.log"
  threads: 1
  conda : "../envs/utility.yml"
  shell:
    """
    rm -rf {params.tmp_dir}/*
    mkdir -p {params.tmp_dir}

    python {params.script_path} \
        --mergedclassify {input.classify_out} \
        --checkvsummary {input.checkv_out} \
        --output {params.tmp_dir}/tmp.csv 2> {log}

    mv {params.tmp_dir}/tmp.csv {output}
    """



rule consensus_filtering:
  name: "viralcontigident.py consensus vOTU filtering"
  input:
    "results/viralcontigident/output/checkv/combined_classification_results.csv"
  output:
    summary = "results/viralcontigident/output/classification_summary_vOTUs.csv", 
    proviruslist = "results/viralcontigident/output/provirus.list.txt",
    viruslist = "results/viralcontigident/output/virus.list.txt"
  params:
    script_path = "workflow/scripts/viralcontigident/consensus_filtering.py",
    genomad_cutoff = config['genomadcutoff'],
    dvf_cutoff = config['dvfcutoff'],
    phamer_pred = config['phamerpred'], 
    tmp_dir = "tmp"
  log: "logs/viralcontigident_consensusfiltering.log"
  threads: 1
  conda: "../envs/utility.yml"
  shell:
    """
    rm -rf {params.tmp_dir}/*
    mkdir -p {params.tmp_dir}

    python {params.script_path} \
        --classification_results {input} \
        --genomad_min_score {params.genomad_cutoff} \
        --dvf_min_score {params.dvf_cutoff} \
        --phamer_pred {params.phamer_pred} \
        --summary_out {params.tmp_dir}/tmp.csv \
        --provirus_list {params.tmp_dir}/tmp.list.1 \
        --virus_list {params.tmp_dir}/tmp.list.2  2> {log}

    mv {params.tmp_dir}/tmp.csv {output.summary}
    mv {params.tmp_dir}/tmp.list.1 {output.proviruslist}
    mv {params.tmp_dir}/tmp.list.2 {output.viruslist}
    """


rule votu:
  name :"viralcontigident.py generate final vOTUs"
  input:
    provirusfasta = "results/viralcontigident/output/checkv/proviruses.fna",
    virusfasta = "results/viralcontigident/output/checkv/viruses.fna", 
    provirushits = "results/viralcontigident/output/provirus.list.txt",
    virushits = "results/viralcontigident/output/virus.list.txt"
  output:
    combined = "results/viralcontigident/output/combined.final.vOTUs.fa",
    provirus = "results/viralcontigident/output/provirus.final.vOTUs.fa",
    virus = "results/viralcontigident/output/virus.final.vOTUs.fa"
  params:
    output_dir = "results/viralcontigident/output/", 
    tmp_dir = "$TMPDIR"
  log: "logs/viralcontigident_vOTUs.log"
  threads: 1
  conda: "../envs/utility.yml"
  shell:
    """
    rm -rf {params.tmp_dir}/*
    mkdir -p {params.tmp_dir} {params.output_dir}

    seqkit replace {input.provirusfasta} --f-use-regexp -p "(.+)_\d\s.+$" -r '$1' | \
        seqkit grep -f {input.provirushits} > {params.tmp_dir}/tmp1.fa 2> {log}
    seqkit replace {input.virusfasta} --f-use-regexp -p "(.+)_\d\s.+$" -r '$1' | \
        seqkit grep -f {input.provirushits} >> {params.tmp_dir}/tmp1.fa 2> {log}

    seqkit grep {input.virusfasta} -f {input.virushits} > {params.tmp_dir}/tmp2.fa 2> {log}
    seqkit grep {input.provirusfasta} -f {input.virushits} >> {params.tmp_dir}/tmp2.fa 2> {log}

    cat {params.tmp_dir}/tmp1.fa {params.tmp_dir}/tmp2.fa > {params.tmp_dir}/tmp3.fa 2> {log}

    mv {params.tmp_dir}/tmp1.fa {output.provirus}
    mv {params.tmp_dir}/tmp2.fa {output.virus}
    mv {params.tmp_dir}/tmp3.fa {output.combined}

    rm -rf {params.tmp_dir}/*

    """

