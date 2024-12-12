import os 
import json 

from rich.console import Console
from rich.progress import Progress
from rich.layout import Layout
from rich.panel import Panel
console = Console()

configdict = config['viral-idenitfy']
logdir=relpath("identify/viral/logs")
benchmarks=relpath("identify/viral/benchmarks")
tmpd=relpath("identify/viral/tmp")

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

n_cores = config['cores'] 
assembler = config['assembler']

### Read single fasta file if input
if config['fasta'] != "":
  fastap = readfasta(config['fasta'])
  sample_id = config["sample-name"]
  assembly_ids = [sample_id]
else:
  samples, assemblies = parse_sample_list(config["samplelist"], datadir, outdir, email, nowstr)
  fastap = relpath(os.path.join("assembly", assembler, "samples/{sample_id}/output/final.contigs.fa"))
  assembly_ids = assemblies.keys()


### MASTER RULE 

rule done_log:
  name: "viral-multitool.smk Done. removing tmp files"
  localrule: True
  input:
    expand(relpath("identify/viral/samples/{sample_id}/intermediate/dvf/final_score.txt"), sample_id=assembly_ids),
    expand(relpath("identify/viral/samples/{sample_id}/intermediate/phamer/out/phamer_prediction.csv"), sample_id=assembly_ids),
    relpath("identify/viral/intermediate/scores/combined.viralcontigs.fa"),
    relpath("identify/viral/intermediate/scores/combined_viral_scores.csv"),
    relpath("identify/viral/output/checkv/combined_classification_results.csv"),
    relpath("identify/viral/output/combined.final.vOTUs.fa"), 
    allhits=relpath("identify/viral/intermediate/scores/combined_allcontigs_scores.csv")
  output:
    os.path.join(logdir, "done.log")
  params:
    filteredcontigs=expand(relpath("identify/viral/samples/{sample_id}/tmp"), sample_id=assembly_ids),
    tmpdir=tmpd
  log: os.path.join(logdir, "done.log")
  shell:
    """
    rm -rf {params.tmpdir}/*
    touch {output}
    """


### RULES

rule filter_contigs:
  name: "viral-multitool.smk filter contigs [length]"
  localrule: True
  input:
    fastap
  output:
    relpath("identify/viral/samples/{sample_id}/tmp/final.contigs.filtered.fa")
  params:
    minlen=configdict['contigminlen'],
    outdir=relpath("identify/viral/samples/{sample_id}/tmp"), 
    tmpdir=os.path.join(tmpd, "contigs/{sample_id}")
  log: os.path.join(logdir, "filtercontig_{sample_id}.log")
  conda: "../envs/utility.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir}/* {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}
    
    seqkit seq {input} --min-len {params.minlen} > {params.tmpdir}/tmp.fa

    mv {params.tmpdir}/tmp.fa {output}
    """



rule genomad_classify:
  name: "viral-multitool.smk geNomad classify" 
  input:
    fna=relpath("identify/viral/samples/{sample_id}/tmp/final.contigs.filtered.fa"),
  output:
    relpath("identify/viral/samples/{sample_id}/intermediate/genomad/final.contigs.filtered_summary/final.contigs.filtered_virus_summary.tsv")
  params:
    genomadparams=configdict['genomadparams'],
    dbdir=configdict['genomaddb'],
    outdir=relpath("identify/viral/samples/{sample_id}/intermediate/genomad/"),
    tmpdir=os.path.join(tmpd, "genomad/{sample_id}")
  log: os.path.join(logdir, "genomad_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "genomad_{sample_id}.log")
  conda: "../envs/genomad.yml"
  threads: min(32, n_cores//3)
  resources:
    mem_mb=lambda wildcards, attempt, input: 24 * 10**3 * attempt
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

    mv {params.tmpdir}/* {params.outdir}
    rm -rf {params.tmpdir}
    """



rule dvf_classify:
  name : "viral-multitool.smk DeepVirFinder classify"
  input:
    fna=relpath("identify/viral/samples/{sample_id}/tmp/final.contigs.filtered.fa")
  output:
    relpath("identify/viral/samples/{sample_id}/intermediate/dvf/final_score.txt")
  params:
    script="workflow/software/DeepVirFinder/dvf.py",
    parameters=configdict['dvfparams'], 
    modeldir="workflow/software/DeepVirFinder/models/",
    outdir=relpath("identify/viral/samples/{sample_id}/intermediate/dvf/"),
    tmpdir=os.path.join(tmpd, "dvf/{sample_id}")
  log: os.path.join(logdir, "dvf_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "dvf_{sample_id}.log")
  conda: "../envs/dvf.yml"
  threads: min(32, n_cores//3)
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: max(6 * threads * 10**3 * attempt, 8000)
  shell:
    """
    rm -rf {params.tmpdir}/* 
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        -i {input.fna} \
        -l 0 \
        -m {params.modeldir} \
        -c {threads} \
        -o {params.tmpdir} \
        {params.parameters} &> {log}

    mv {params.tmpdir}/* {output}
    rm -rf {params.tmpdir} 
    """



rule phamer_classify:
  name: "viral-multitool.smk PhaMer classify"
  input:
    fna=relpath("identify/viral/samples/{sample_id}/tmp/final.contigs.filtered.fa")
  output:
    relpath("identify/viral/samples/{sample_id}/intermediate/phamer/final_prediction/phamer_prediction.tsv")
  params:
    parameters=configdict['phamerparams'],
    dbdir=configdict['phamerdb'],
    outdir=relpath("identify/viral/samples/{sample_id}/intermediate/phamer/"),
    tmpdir=os.path.join(tmpd, "phamer/{sample_id}")
  log: os.path.join(logdir, "phamer_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "phamer_{sample_id}.log")
  conda: "../envs/phabox2.yml"
  threads: min(32, n_cores//3)
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: max(1 * threads * 10**3 * attempt, 8000)
  shell:
    """
    rm -rf {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    phabox2 --task phamer \
        --contigs {input.fna} \
        --len 0 \
        --threads {threads} \
        --outpth {params.tmpdir} \
        --dbdir {params.dbdir} \
        {params.parameters} &> {log}

    mv -f {params.tmpdir}/* {params.outdir}
    rm -rf {params.tmpdir}
    """


# NOTE, THE PYTHON CODE NEEDS TO BE CHANGED AS PHAMERS OUTPUT IS NOW A TSV ONE

rule merge_outputs:
  name : "viral-multitool.smk merge classification outputs"
  input:
    genomadout=relpath("identify/viral/samples/{sample_id}/intermediate/genomad/final.contigs.filtered_summary/final.contigs.filtered_virus_summary.tsv"),
    dvfout=relpath("identify/viral/samples/{sample_id}/intermediate/dvf/final_score.txt"), 
    phamerout=relpath("identify/viral/samples/{sample_id}/intermediate/phamer/out/phamer_prediction.tsv"),
  output:
    relpath("identify/viral/samples/{sample_id}/output/merged_scores.csv")
  params:
    script="workflow/scripts/identify/viral/mergeout.py",
    genomadminlen=configdict['genomadminlen'],
    dvfminlen=configdict['dvfminlen'], 
    phamerminlen=configdict['phamerminlen'], 
    outdir=relpath("identify/viral/samples/{sample_id}/output/"),
    tmpdir=os.path.join(tmpd, "merge/{sample_id}")
  log: os.path.join(logdir, "mergeout_{sample_id}.log")
  conda: "../envs/seqkit-biopython.yml"
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
  name: "viral-multitool.smk filter viral contigs"
  input:
    fna=relpath("identify/viral/samples/{sample_id}/tmp/final.contigs.filtered.fa"),
    scrs=relpath("identify/viral/samples/{sample_id}/output/merged_scores.csv")
  output:
    fna=relpath("identify/viral/samples/{sample_id}/output/viral.contigs.fa"),
    scrs=relpath("identify/viral/samples/{sample_id}/output/merged_scores_filtered.csv"),
    hits=relpath("identify/viral/samples/{sample_id}/output/viralhits_list")
  params:
    script="workflow/scripts/identify/viral/filtercontig_scores.py",
    genomad_cutoff=configdict['genomadcutoff_p'], 
    dvf_cutoff=configdict['dvfcutoff_p'], 
    dvf_pvalmax=configdict['dvfpval_p'],
    phamer_cutoff=configdict['phamercutoff_p'], 
    phamer_pred=configdict['phamerpred_p'], 
    tmpdir=os.path.join(tmpd, "filter/{sample_id}")
  log: os.path.join(logdir, "filteroutput_{sample_id}.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}

    python {params.script} \
        --csv_path {input.scrs} \
        --genomad_min_score {params.genomad_cutoff} \
        --dvf_min_score {params.dvf_cutoff} \
        --dvf_max_pval {params.dvf_pvalmax} \
        --phamer_pred {params.phamer_pred} \
        --phamer_min_score {params.phamer_cutoff} \
        --output_path {params.tmpdir}/tmp.csv \
        --hitlist_path {params.tmpdir}/tmplist

    cat {params.tmpdir}/tmplist | uniq > {output.hits}
    mv {params.tmpdir}/tmp.csv {output.scrs}
    
    seqkit grep {input.fna} -f {output.hits} | seqkit replace -p  "\s.*" -r "" | seqkit replace -p $ -r _{wildcards.sample_id}  > {params.tmpdir}/tmp.fa 2> {log}
    mv {params.tmpdir}/tmp.fa {output.fna}

    rm -r {params.tmpdir}
    """
 

rule cat_contigs:
  name : "viral-multitool.smk combine viral contigs"
  input:
    fna=expand(relpath("identify/viral/samples/{sample_id}/output/viral.contigs.fa"), sample_id=assembly_ids),
    scrs=expand(relpath("identify/viral/samples/{sample_id}/output/merged_scores_filtered.csv"), sample_id=assembly_ids),
    allhits=expand(relpath("identify/viral/samples/{sample_id}/output/merged_scores.csv"), sample_id=assembly_ids)
  output: 
    fna=relpath("identify/viral/intermediate/scores/combined.viralcontigs.fa"),
    scrs=relpath("identify/viral/intermediate/scores/combined_viral_scores.csv"),
    allhits=relpath("identify/viral/intermediate/scores/combined_allcontigs_scores.csv")
  params:
    script="workflow/scripts/identify/viral/mergeout_scores.py", 
    names=list(assembly_ids),
    tmpdir=tmpd
  log: os.path.join(logdir, "catcontigs.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  resources:
    maxcores=1
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}

    echo "{params.names}" > {params.tmpdir}/tmp.names
    echo "{input.scrs}" > {params.tmpdir}/tmp.hits.paths
    echo "{input.allhits}" > {params.tmpdir}/tmp.allhits.paths

    python {params.script} \
        --names {params.tmpdir}/tmp.names \
        --csvs {params.tmpdir}/tmp.hits.paths > {params.tmpdir}/tmp.csv 2> {log}
        
    python {params.script} \
        --names {params.tmpdir}/tmp.names \
        --csvs {params.tmpdir}/tmp.allhits.paths > {params.tmpdir}/tmp.allhits.csv 2> {log}

    cat {input.fna} > {params.tmpdir}/tmp.fa 2> {log}

    mv {params.tmpdir}/tmp.csv {output.scrs}
    mv {params.tmpdir}/tmp.allhits.csv {output.allhits}
    mv {params.tmpdir}/tmp.fa {output.fna}

    rm -rf {params.tmpdir}
    """


# THEN GOES INTO 
# 1) VIRAL BINNING [optional]
# 2) CLUSTERING [sensitive or fast]
# 3) CHECKV-PYHMMER 
# THEN COMES BACK HERE TO GET vCONTIGS


rule combine_classifications:
  name: "viral-multitool.smk combine derepped classification results"
  input:
    checkv_out=relpath("identify/viral/output/checkv/quality_summary.tsv"),
    classify_out=relpath("identify/viral/intermediate/scores/combined_viral_scores.csv")
  output:
    relpath("identify/viral/output/checkv/combined_classification_results.csv")
  params:
    script="workflow/scripts/identify/viral/combineclassify.py",
    tmpdir=tmpd
  log: os.path.join(logdir, "combine_classification.log")
  threads: 1
  conda : "../envs/seqkit-biopython.yml"
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
  name: "viral-multitool.smk consensus vOTU filtering"
  input:
    relpath("identify/viral/output/checkv/combined_classification_results.csv")
  output:
    summary=relpath("identify/viral/output/classification_summary_vOTUs.csv"),
    proviruslist=relpath("identify/viral/output/provirus.list.txt"),
    viruslist=relpath("identify/viral/output/virus.list.txt")
  params:
    script="workflow/scripts/identify/viral/consensus_filtering.py",
    genomad=configdict['genomadcutoff_s'],
    dvf=configdict['dvfcutoff_s'],
    phamer=configdict['phamercutoff_s'], 
    tmpdir=tmpd
  log: os.path.join(logdir, "consensus_filtering.log")
  threads: 1
  conda: "../envs/seqkit-biopython.yml"
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir}

    python {params.script} \
        --classification_results {input} \
        --genomad_min_score {params.genomad} \
        --dvf_min_score {params.dvf} \
        --phamer_min_score {params.phamer} \
        --summary_out {params.tmpdir}/tmp.csv \
        --provirus_list {params.tmpdir}/tmp.list.1 \
        --virus_list {params.tmpdir}/tmp.list.2  2> {log}

    mv {params.tmpdir}/tmp.csv {output.summary}
    mv {params.tmpdir}/tmp.list.1 {output.proviruslist}
    mv {params.tmpdir}/tmp.list.2 {output.viruslist}
    """



rule votu:
  name: "viral-multitool.smk generate final vOTUs"
  input:
    provirusfasta=relpath("identify/viral/output/checkv/proviruses.fna"),
    virusfasta=relpath("identify/viral/output/checkv/viruses.fna"), 
    provirushits=relpath("identify/viral/output/provirus.list.txt"),
    virushits=relpath("identify/viral/output/virus.list.txt")
  output:
    combined=relpath("identify/viral/output/combined.final.vOTUs.fa"),
    provirus=relpath("identify/viral/output/provirus.final.vOTUs.fa"),
    virus=relpath("identify/viral/output/virus.final.vOTUs.fa"), 
    tsv=relpath("identify/viral/output/GC_content_vOTUs.tsv")
  params:
    outdir=relpath("identify/viral/output/"),
    tmpdir=tmpd
  log: os.path.join(logdir, "vOTUs.log")
  threads: 1
  conda: "../envs/seqkit-biopython.yml"
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir} {params.outdir}

    seqkit replace {input.provirusfasta} --f-use-regexp -p "(.+)_\d\s.+$" -r '$1' \
        | seqkit grep -f {input.provirushits} > {params.tmpdir}/tmp1.fa 2> {log}
    seqkit replace {input.virusfasta} --f-use-regexp -p "(.+)_\d\s.+$" -r '$1' \
        | seqkit grep -f {input.provirushits} >> {params.tmpdir}/tmp1.fa 2> {log}

    seqkit grep {input.virusfasta} -f {input.virushits} > {params.tmpdir}/tmp2.fa 2> {log}
    seqkit grep {input.provirusfasta} -f {input.virushits} >> {params.tmpdir}/tmp2.fa 2> {log}

    cat {params.tmpdir}/tmp1.fa {params.tmpdir}/tmp2.fa  > {params.tmpdir}/tmp3.fa 2> {log}

    seqkit rmdup {params.tmpdir}/tmp1.fa > {output.provirus} 2> {log}
    seqkit rmdup {params.tmpdir}/tmp2.fa > {output.virus} 2> {log}
    seqkit rmdup {params.tmpdir}/tmp3.fa > {output.combined} 2> {log}
    
    echo -e "seq_id\tlength\tGC_percent" > {params.tmpdir}/tmp.tsv
    seqkit fx2tab -g -l -n {output.combined} >> {params.tmpdir}/tmp.tsv
    mv {params.tmpdir}/tmp.tsv {output.tsv}

    rm -rf {params.tmpdir}/*
    """

