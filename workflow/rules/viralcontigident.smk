configfile: "config/viralcontigident.yml"


rule genomad_classify:
  name: "viralcontigident.py geNomad classify" 
  input:
    fna=os.path.join(config['contigdir'], "{sample_id}/output/final.contigs.fa"),
  output:
    "results/viralcontigident/samples/{sample_id}/intermediate/genomad/final.contigs_summary/final.contigs_virus_summary.tsv"
  params:
    genomadparams=config['genomadparams'],
    dbdir=config['genomaddb'],
    outdir="results/viralcontigident/samples/{sample_id}/intermediate/genomad/",
    tmpdir="$TMPDIR/{sample_id}"
  log: "logs/viralcontigident_{sample_id}_genomad.log"
  benchmark: "benchmarks/viralcontigident_{sample_id}_genomad.log"
  conda: "../envs/genomad.yml"
  threads: 8
  resources:
    #runtime_min=lambda wildcards, attempt: attempt*attempt*120
    mem_mb=lambda wildcards, attempt: attempt * 72 * 10**3
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir} 2> {log}
    mkdir -p {params.tmpdir} {params.outdir} 2> {log}

    genomad end-to-end \
        {input.fna} \
        {params.tmpdir} \
        {params.dbdir} \
        --threads {threads} \
        --cleanup \
        {params.genomadparams} &> {log}

    mv {params.tmpdir}/* {params.outdir} 2> {log}
    rm -rf {params.tmpdir} 2> {log}
    """



rule dvf_classify:
  name : "viralcontigident.py DeepVirFinder classify"
  input:
    fna=os.path.join(config['contigdir'], "{sample_id}/output/final.contigs.fa")
  output:
    "results/viralcontigident/samples/{sample_id}/intermediate/dvf/final_score.txt"
  params:
    script="workflow/software/DeepVirFinder/dvf.py",
    parameters=config['dvfparams'], 
    modeldir="workflow/software/DeepVirFinder/models/",
    outdir="results/viralcontigident/samples/{sample_id}/intermediate/dvf/",
    tmpdir="$TMPDIR/{sample_id}"
  log: "logs/viralcontigident_{sample_id}_dvf.log"
  benchmark: "benchmarks/viralcontigident_{sample_id}_dvf.log"
  conda: "../envs/dvf.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
  shell:
    """
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        -i {input.fna} \
        -m {params.modeldir} \
        -c {threads} \
        -o {params.tmpdir} \
        {params.parameters} &> {log}

    mv {params.tmpdir}/* {output}
    rm -rf {params.tmpdir}
    """


rule phamer_classify:
  name: "viralcontigident.py PhaMer classify"
  input:
    fna=os.path.join(config['contigdir'], "{sample_id}/output/final.contigs.fa")
  output:
    "results/viralcontigident/samples/{sample_id}/intermediate/phamer/out/phamer_prediction.csv"
  params:
    script="workflow/software/PhaBOX/PhaMer_single.py",
    scriptdir="workflow/software/PhaBOX/scripts/",
    parameters=config['phamerparams'],
    paramsdir="workflow/params/phabox/",
    dbdir=config['phamerdb'],
    outdir="results/viralcontigident/samples/{sample_id}/intermediate/phamer/",
    tmpdir="$TMPDIR/{sample_id}"
  log: "logs/viralcontigident_{sample_id}_phamer.log"
  benchmark: "benchmarks/viralcontigident_{sample_id}_phamer.log"
  conda: "../envs/phabox.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
  shell:
    """
    rm -rf {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        --contigs {input.fna} \
        --threads {threads} \
        --rootpth {params.tmpdir} \
        --dbdir {params.dbdir} \
        --parampth {params.paramsdir} \
        --scriptpth {params.scriptdir} \
        {params.parameters} &> {log}

    mv -f {params.tmpdir}/* {params.outdir}
    rm -rf {params.tmpdir}
    """


rule merge_outputs:
  name : "viralcontigident.py merge classification outputs"
  input:
    genomadout="results/viralcontigident/samples/{sample_id}/intermediate/genomad/final.contigs_summary/final.contigs_virus_summary.tsv",
    dvfout="results/viralcontigident/samples/{sample_id}/intermediate/dvf/final_score.txt", 
    phamerout="results/viralcontigident/samples/{sample_id}/intermediate/phamer/out/phamer_prediction.csv",
  output:
    "results/viralcontigident/samples/{sample_id}/output/merged_scores.csv"
  params:
    script="workflow/scripts/viralcontigident/mergeout.py",
    genomadminlen=config['genomadminlen'],
    dvfminlen=config['dvfminlen'], 
    phamerminlen=config['dvfminlen'], 
    outdir="results/viralcontigident/samples/{sample_id}/output/",
    tmpdir="$TMPDIR/{sample_id}"
  log: "logs/viralcontigident_{sample_id}_mergeoutput.log"
  conda: "../envs/utility.yml"
  threads: 1
  shell:
    """
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        --genomadout {input.genomadout} \
        --dvfout {input.dvfout} \
        --phamerout {input.phamerout} \
        --genomadminlen {params.genomadminlen} \
        --dvfminlen {params.dvfminlen} \
        --phamerminlen {params.phamerminlen} \
        --output {params.tmpdir}/tmp.csv &> {log}

    mv {params.tmpdir}/* {output}
    rm -rf {params.tmpdir}
    """
    

rule filter_outputs:
  name : "viralcontigident.py filter viral contigs"
  input:
    contig_file=os.path.join(config['contigdir'], "{sample_id}/output/final.contigs.fa"),
    merged_scrs="results/viralcontigident/samples/{sample_id}/output/merged_scores.csv"
  output:
    filtered_contigs="results/viralcontigident/samples/{sample_id}/output/viral.contigs.fa",
    filtered_scrs="results/viralcontigident/samples/{sample_id}/output/merged_scores_filtered.csv",
    positive_hits="results/viralcontigident/samples/{sample_id}/output/viralhits_list"
  params:
    script="workflow/scripts/viralcontigident/filtercontig_scores.py",
    genomad_cutoff=config['genomadcutoff'], 
    dvf_cutoff=config['dvfcutoff'], 
    dvf_pvalmax=config['dvfpval'],
    phamer_cutoff=config['phamercutoff'], 
    phamer_pred=config['phamerpred'], 
    tmpdir="$TMPDIR/{sample_id}"
  log: "logs/viralcontigident_{sample_id}_filtercontigs.log"
  conda: "../envs/utility.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}

    python {params.script} \
        --csv_path {input.merged_scrs} \
        --genomad_min_score {params.genomad_cutoff} \
        --dvf_min_score {params.dvf_cutoff} \
        --dvf_max_pval {params.dvf_pvalmax} \
        --phamer_pred {params.phamer_pred} \
        --phamer_min_score {params.phamer_cutoff} \
        --output_path {params.tmpdir}/tmp.csv \
        --hitlist_path {params.tmpdir}/tmplist

    cat {params.tmpdir}/tmplist | uniq > {output.positive_hits}
    mv {params.tmpdir}/tmp.csv {output.filtered_scrs}
    
    seqkit grep {input.contig_file} -f {output.positive_hits} | seqkit replace -p  "\s.*" -r "" |seqkit replace -p $ -r _{wildcards.sample_id}  > {params.tmpdir}/tmp.fa 2> {log}
    mv {params.tmpdir}/tmp.fa {output.filtered_contigs}

    rm -r {params.tmpdir}
    """
 

rule cat_contigs:
  name : "viralcontigident.py combine viral contigs"
  input:
    fasta=expand("results/viralcontigident/samples/{sample_id}/output/viral.contigs.fa", sample_id = samples.keys()),
    scores=expand("results/viralcontigident/samples/{sample_id}/output/merged_scores_filtered.csv", sample_id = samples.keys())
  output: 
    fasta="results/viralcontigident/intermediate/scores/combined.viralcontigs.fa",
    scores="results/viralcontigident/intermediate/scores/combined_viral_scores.csv"
  params:
    script="workflow/scripts/viralcontigident/mergeout_scores.py", 
    tmpdir="$TMPDIR/viralcontigident"
  log: "logs/viralcontigident_catcontigs.log"
  conda: "../envs/utility.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}

    python {params.script} {input.scores} > {params.tmpdir}/tmp.csv 2> {log}
    cat {input.fasta} > {params.tmpdir}/tmp.fa 2> {log}

    mv {params.tmpdir}/tmp.csv {output.scores}
    mv {params.tmpdir}/tmp.fa {output.fasta}

    rm -rf {params.tmpdir}
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
    outdir="results/viralcontigident/output/checkv",
    tmpdir="$TMPDIR/checkv",
    dbdir="workflow/database/checkv"
  log: "logs/viralcontigident_checkv.log"
  benchmark: "benchmarks/viralcontigident_checkv.log"
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 72 * 10**3
  conda: "../envs/checkv.yml"
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}/*
    mkdir -p {params.tmpdir} {params.outdir}

    checkv end_to_end {input} {params.tmpdir} -d {params.dbdir} -t {threads} {params.checkvparams} 2> {log}

    mv {params.tmpdir}/* {params.outdir}
    rm -rf {params.tmpdir}"""


rule combine_classifications:
  name: "viralcontigident.py combine derepped classification results"
  input:
    checkv_out="results/viralcontigident/output/checkv/quality_summary.tsv",
    classify_out="results/viralcontigident/intermediate/scores/combined_viral_scores.csv"
  output:
    "results/viralcontigident/output/checkv/combined_classification_results.csv"
  params:
    script="workflow/scripts/viralcontigident/combineclassify.py",
    tmpdir="$TMPDIR"
  log: "logs/viralcontigident_combineclassification.log"
  threads: 1
  conda : "../envs/utility.yml"
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir}

    python {params.script} \
        --mergedclassify {input.classify_out} \
        --checkvsummary {input.checkv_out} \
        --output {params.tmpdir}/tmp.csv 2> {log}

    mv {params.tmpdir}/tmp.csv {output}
    """



rule consensus_filtering:
  name: "viralcontigident.py consensus vOTU filtering"
  input:
    "results/viralcontigident/output/checkv/combined_classification_results.csv"
  output:
    summary="results/viralcontigident/output/classification_summary_vOTUs.csv", 
    proviruslist="results/viralcontigident/output/provirus.list.txt",
    viruslist="results/viralcontigident/output/virus.list.txt"
  params:
    script="workflow/scripts/viralcontigident/consensus_filtering.py",
    genomad=config['genomadcutoff'],
    dvf=config['dvfcutoff'],
    phamer=config['phamerpred'], 
    tmpdir="tmp"
  log: "logs/viralcontigident_consensusfiltering.log"
  threads: 1
  conda: "../envs/utility.yml"
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir}

    python {params.script} \
        --classification_results {input} \
        --genomad_min_score {params.genomad} \
        --dvf_min_score {params.dvf} \
        --phamer_pred {params.phamer} \
        --summary_out {params.tmpdir}/tmp.csv \
        --provirus_list {params.tmpdir}/tmp.list.1 \
        --virus_list {params.tmpdir}/tmp.list.2  2> {log}

    mv {params.tmpdir}/tmp.csv {output.summary}
    mv {params.tmpdir}/tmp.list.1 {output.proviruslist}
    mv {params.tmpdir}/tmp.list.2 {output.viruslist}
    """


rule votu:
  name :"viralcontigident.py generate final vOTUs"
  input:
    provirusfasta="results/viralcontigident/output/checkv/proviruses.fna",
    virusfasta="results/viralcontigident/output/checkv/viruses.fna", 
    provirushits="results/viralcontigident/output/provirus.list.txt",
    virushits="results/viralcontigident/output/virus.list.txt"
  output:
    combined="results/viralcontigident/output/combined.final.vOTUs.fa",
    provirus="results/viralcontigident/output/provirus.final.vOTUs.fa",
    virus="results/viralcontigident/output/virus.final.vOTUs.fa"
  params:
    outdir="results/viralcontigident/output/", 
    tmpdir="$TMPDIR"
  log: "logs/viralcontigident_vOTUs.log"
  threads: 1
  conda: "../envs/utility.yml"
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir} {params.outdir}

    seqkit replace {input.provirusfasta} --f-use-regexp -p "(.+)_\d\s.+$" -r '$1' | \
        seqkit grep -f {input.provirushits} > {params.tmpdir}/tmp1.fa 2> {log}
    seqkit replace {input.virusfasta} --f-use-regexp -p "(.+)_\d\s.+$" -r '$1' | \
        seqkit grep -f {input.provirushits} >> {params.tmpdir}/tmp1.fa 2> {log}

    seqkit grep {input.virusfasta} -f {input.virushits} > {params.tmpdir}/tmp2.fa 2> {log}
    seqkit grep {input.provirusfasta} -f {input.virushits} >> {params.tmpdir}/tmp2.fa 2> {log}

    cat {params.tmpdir}/tmp1.fa {params.tmpdir}/tmp2.fa > {params.tmpdir}/tmp3.fa 2> {log}

    mv {params.tmpdir}/tmp1.fa {output.provirus}
    mv {params.tmpdir}/tmp2.fa {output.virus}
    mv {params.tmpdir}/tmp3.fa {output.combined}

    rm -rf {params.tmpdir}/*
    """

