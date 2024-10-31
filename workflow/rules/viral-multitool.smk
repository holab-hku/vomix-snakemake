import os 
import json 

from rich.console import Console
from rich.progress import Progress
from rich.layout import Layout
from rich.panel import Panel
console = Console()


configdict = config['viral-contigident']
logdir=relpath("viralcontigident/logs")
benchmarks=relpath("viralcontigident/benchmarks")
tmpd=relpath("viralcontigident/tmp")

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)


if isinstance(config['cores'], int):
 n_cores = config['cores'] 
else:
  console.print(Panel.fit(f"config['cores'] is not an integer: {config['cores']}, you can change the parameter in config/config.yml file", title="Error", subtitle="config['cores'] not integer"))
  sys.exit(1)


############################
# Single-Sample Processing #
############################

# the "contigfile" is not the nested one

if config['inputdir']!="":

  indir = cleanpath(config['inputdir'])
  console.print(f"\nconfig['inputdir'] not empty, using '{indir}' as input for viral contig identification.")
  console.print("File names without the .fa extension will be used as sample IDs.")
  cwd = os.getcwd()
  indir_path = os.path.join(cwd, indir)

  if not os.path.exists(indir):
    console.print(Panel.fit(f"The input file path '{indir}' does not exist.", title="Error", subtitle="Contig Directory"))
    sys.exit(1)

  fasta_files = [f for f in os.listdir(indir_path) if f.endswith('.fa')]
  if len(fasta_files) == 0:
    console.print(Panel.fit(f"There are no files ending with .fa in '{indir_path}', other fasta extensions are not accepted (for now) :(", title="Error", subtitle="No .fa Files"))
    sys.exit(1) 

  assembly_ids = [os.path.basename(fasta_file).rsplit(".", 1)[0] for fasta_file in fasta_files]
  wildcards_p = os.path.join(indir, "{sample_id}.fa") 
  outdir_p = os.path.join(cwd, relpath("viralcontigident/"))
  console.print(f"Creating output directory: '{outdir_p}'.\n")

else:
  wildcards_p = relpath("assembly/samples/{sample_id}/output/final.contigs.fa")
  assembly_ids = assemblies.keys()



###########################
# Multi-sample Processing #
###########################


### MASTER RULE 

rule done_log:
  name: "viral-multitool.smk Done. removing tmp files"
  localrule: True
  input:
    expand(relpath("viralcontigident/samples/{sample_id}/intermediate/dvf/final_score.txt"), sample_id=assembly_ids),
    expand(relpath("viralcontigident/samples/{sample_id}/intermediate/phamer/out/phamer_prediction.csv"), sample_id=assembly_ids),
    relpath("viralcontigident/intermediate/scores/combined.viralcontigs.fa"),
    relpath("viralcontigident/intermediate/scores/combined_viral_scores.csv"),
    relpath("viralcontigident/output/checkv/combined_classification_results.csv"),
    relpath("viralcontigident/output/combined.final.vOTUs.fa"), 
    allhits=relpath("viralcontigident/intermediate/scores/combined_allcontigs_scores.csv")
  output:
    os.path.join(logdir, "done.log")
  params:
    filteredcontigs=expand(relpath("viralcontigident/samples/{sample_id}/tmp"), sample_id=assembly_ids),
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
    wildcards_p
  output:
    relpath("viralcontigident/samples/{sample_id}/tmp/final.contigs.filtered.fa")
  params:
    minlen=configdict['contigminlen'],
    outdir=relpath("viralcontigident/samples/{sample_id}/tmp"), 
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
    fna=relpath("viralcontigident/samples/{sample_id}/tmp/final.contigs.filtered.fa"),
  output:
    relpath("viralcontigident/samples/{sample_id}/intermediate/genomad/final.contigs.filtered_summary/final.contigs.filtered_virus_summary.tsv")
  params:
    genomadparams=configdict['genomadparams'],
    dbdir=configdict['genomaddb'],
    outdir=relpath("viralcontigident/samples/{sample_id}/intermediate/genomad/"),
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
    fna=relpath("viralcontigident/samples/{sample_id}/tmp/final.contigs.filtered.fa")
  output:
    relpath("viralcontigident/samples/{sample_id}/intermediate/dvf/final_score.txt")
  params:
    script="workflow/software/DeepVirFinder/dvf.py",
    parameters=configdict['dvfparams'], 
    modeldir="workflow/software/DeepVirFinder/models/",
    outdir=relpath("viralcontigident/samples/{sample_id}/intermediate/dvf/"),
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
    fna=relpath("viralcontigident/samples/{sample_id}/tmp/final.contigs.filtered.fa")
  output:
    relpath("viralcontigident/samples/{sample_id}/intermediate/phamer/out/phamer_prediction.csv")
  params:
    script="workflow/software/PhaBOX/PhaMer_single.py",
    scriptdir="workflow/software/PhaBOX/scripts/",
    parameters=configdict['phamerparams'],
    paramsdir="workflow/params/phabox/",
    dbdir=configdict['phamerdb'],
    outdir=relpath("viralcontigident/samples/{sample_id}/intermediate/phamer/"),
    tmpdir=os.path.join(tmpd, "phamer/{sample_id}")
  log: os.path.join(logdir, "phamer_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "phamer_{sample_id}.log")
  conda: "../envs/phabox.yml"
  threads: min(32, n_cores//3)
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: max(1 * threads * 10**3 * attempt, 8000)
  shell:
    """
    rm -rf {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        --contigs {input.fna} \
        --len 0 \
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
  name : "viral-multitool.smk merge classification outputs"
  input:
    genomadout=relpath("viralcontigident/samples/{sample_id}/intermediate/genomad/final.contigs.filtered_summary/final.contigs.filtered_virus_summary.tsv"),
    dvfout=relpath("viralcontigident/samples/{sample_id}/intermediate/dvf/final_score.txt"), 
    phamerout=relpath("viralcontigident/samples/{sample_id}/intermediate/phamer/out/phamer_prediction.csv"),
  output:
    relpath("viralcontigident/samples/{sample_id}/output/merged_scores.csv")
  params:
    script="workflow/scripts/viralcontigident/mergeout.py",
    genomadminlen=configdict['genomadminlen'],
    dvfminlen=configdict['dvfminlen'], 
    phamerminlen=configdict['dvfminlen'], 
    outdir=relpath("viralcontigident/samples/{sample_id}/output/"),
    tmpdir=os.path.join(tmpd, "merge/{sample_id}")
  log: os.path.join(logdir, "mergeout_{sample_id}.log")
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
  name: "viral-multitool.smk filter viral contigs"
  input:
    fna=relpath("viralcontigident/samples/{sample_id}/tmp/final.contigs.filtered.fa"),
    scrs=relpath("viralcontigident/samples/{sample_id}/output/merged_scores.csv")
  output:
    fna=relpath("viralcontigident/samples/{sample_id}/output/viral.contigs.fa"),
    scrs=relpath("viralcontigident/samples/{sample_id}/output/merged_scores_filtered.csv"),
    hits=relpath("viralcontigident/samples/{sample_id}/output/viralhits_list")
  params:
    script="workflow/scripts/viralcontigident/filtercontig_scores.py",
    genomad_cutoff=configdict['genomadcutoff_p'], 
    dvf_cutoff=configdict['dvfcutoff_p'], 
    dvf_pvalmax=configdict['dvfpval_p'],
    phamer_cutoff=configdict['phamercutoff_p'], 
    phamer_pred=configdict['phamerpred_p'], 
    tmpdir=os.path.join(tmpd, "filter/{sample_id}")
  log: os.path.join(logdir, "filteroutput_{sample_id}.log")
  conda: "../envs/utility.yml"
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
    fna=expand(relpath("viralcontigident/samples/{sample_id}/output/viral.contigs.fa"), sample_id=assembly_ids),
    scrs=expand(relpath("viralcontigident/samples/{sample_id}/output/merged_scores_filtered.csv"), sample_id=assembly_ids),
    allhits=expand(relpath("viralcontigident/samples/{sample_id}/output/merged_scores.csv"), sample_id=assembly_ids)
  output: 
    fna=relpath("viralcontigident/intermediate/scores/combined.viralcontigs.fa"),
    scrs=relpath("viralcontigident/intermediate/scores/combined_viral_scores.csv"),
    allhits=relpath("viralcontigident/intermediate/scores/combined_allcontigs_scores.csv")
  params:
    script="workflow/scripts/viralcontigident/mergeout_scores.py", 
    names=list(assembly_ids),
    tmpdir=tmpd
  log: os.path.join(logdir, "catcontigs.log")
  conda: "../envs/utility.yml"
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
    checkv_out=relpath("viralcontigident/output/checkv/quality_summary.tsv"),
    classify_out=relpath("viralcontigident/intermediate/scores/combined_viral_scores.csv")
  output:
    relpath("viralcontigident/output/checkv/combined_classification_results.csv")
  params:
    script="workflow/scripts/viralcontigident/combineclassify.py",
    tmpdir=tmpd
  log: os.path.join(logdir, "combine_classification.log")
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
  name: "viral-multitool.smk consensus vOTU filtering"
  input:
    relpath("viralcontigident/output/checkv/combined_classification_results.csv")
  output:
    summary=relpath("viralcontigident/output/classification_summary_vOTUs.csv"),
    proviruslist=relpath("viralcontigident/output/provirus.list.txt"),
    viruslist=relpath("viralcontigident/output/virus.list.txt")
  params:
    script="workflow/scripts/viralcontigident/consensus_filtering.py",
    genomad=configdict['genomadcutoff_s'],
    dvf=configdict['dvfcutoff_s'],
    phamer=configdict['phamercutoff_s'], 
    tmpdir=tmpd
  log: os.path.join(logdir, "consensus_filtering.log")
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
    provirusfasta=relpath("viralcontigident/output/checkv/proviruses.fna"),
    virusfasta=relpath("viralcontigident/output/checkv/viruses.fna"), 
    provirushits=relpath("viralcontigident/output/provirus.list.txt"),
    virushits=relpath("viralcontigident/output/virus.list.txt")
  output:
    combined=relpath("viralcontigident/output/combined.final.vOTUs.fa"),
    provirus=relpath("viralcontigident/output/provirus.final.vOTUs.fa"),
    virus=relpath("viralcontigident/output/virus.final.vOTUs.fa"), 
    tsv=relpath("viralcontigident/output/GC_content_vOTUs.tsv")
  params:
    outdir=relpath("viralcontigident/output/"),
    tmpdir=tmpd
  log: os.path.join(logdir, "vOTUs.log")
  threads: 1
  conda: "../envs/utility.yml"
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

