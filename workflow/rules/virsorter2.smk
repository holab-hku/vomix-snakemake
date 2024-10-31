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



############################
# Single-Sample Processing #
############################

if config["contigfile"] != "":
  path=config["contigfile"]
  cwd=os.getcwd()
  path_full=os.path.join(cwd, path)

  if not os.path.exists(path):
    console.print(Panel.fit("The contig file path provided does not exist.", title="Error", subtitle="Contig File Path"))
    sys.exit(1)

  assembly_id=os.path.basename(path).rsplit(".", 1)[0]
  outdir_p=os.path.join(cwd, relpath("assembly/samples"), assembly_id, "output/")
  outfile=os.path.join(outdir_p, "final.contigs.fa")

  os.makedirs(outdir_p, exist_ok=True)

  try:
    os.symlink(path_full, outfile)
  except FileExistsError:
    console.print(Panel.fit(f"Warning:[dim] Contig path already exists '{outdir}final.contigs.fa'. The existing file will be used for viral contig annotation. If needed, please delete the current contig file path and try again.", title="Warning", subtitle="Contig Path Exists"))
    pass

  assembly_ids=[assembly_id]

else:
  assembly_ids=assemblies.keys()



###########################
# Multi-sample Processing #
###########################


### MASTER RULE 

rule done_log:
  name: "viral-contigident.py Done. removing tmp files"
  localrule: True
  input:
    expand(relpath("viralcontigident/samples/{sample_id}/intermediate/dvf/final_score.txt"), sample_id=assembly_ids),
    expand(relpath("viralcontigident/samples/{sample_id}/intermediate/phamer/out/phamer_prediction.csv"), sample_id=assembly_ids),
    expand(relpath("viralcontigident/samples/{sample_id}/intermediate/vs2/final-viral-score.tsv"), sample_id=assembly_ids),
    relpath("viralcontigident/intermediate/scores/combined.viralcontigs.fa"),
    relpath("viralcontigident/intermediate/scores/combined_viral_scores.csv"),
    relpath("viralcontigident/output/checkv/combined_classification_results.csv"),
    relpath("viralcontigident/output/combined.final.vOTUs.fa")
  output:
    os.path.join(logdir, "done.log")
  params:
    filteredcontigs=expand(relpath("viralcontigident/samples/{sample_id}/tmp"), sample_id=assembly_ids),
    tmpdir=tmpd
  log: os.path.join(logdir, "done.log")
  shell:
    """
    rm -rf {params.filteredcontigs} {params.tmpdir}/*
    touch {output}
    """


### RULES

rule filter_contigs:
  name: "viral-contigident.py filter contigs [length]"
  localrule: True
  input:
    relpath("assembly/samples/{sample_id}/output/final.contigs.fa")
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


rule vs2_classify:
  name: "viral-contigident.py virsorter2 classify (--virsorter2)"
  input:
    fna=relpath("viralcontigident/samples/{sample_id}/tmp/final.contigs.filtered.fa"),
  output:
    relpath("viralcontigident/samples/{sample_id}/intermediate/vs2/final-viral-score.tsv")
  params:
    vs2params=configdict['vs2params'],
    dbdir=configdict['vs2db'],
    outdir=relpath("viralcontigident/samples/{sample_id}/intermediate/vs2/"), 
    tmpdir=os.path.join(tmpd, "vs2/{sample_id}")
  log: os.path.join(logdir, "vs2_{sample_id}.log")
  benchmark: os.path.join(logdir, "vs2_{sample_id}.log")
  conda: "../envs/vs2.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt, input: 24 * 10**3 * attempt + (input.size_mb * 10)
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    virsorter run \
        -i {input.fna} \
        -w {params.tmpdir} \
        --db-dir {params.dbdir} \
       -j {threads} \
        {params.vs2params} \
        all &> {log}
  
    mv {params.tmpdir}/* {params.output_dir}
  """


rule dvf_classify:
  name : "viral-contigident.py DeepVirFinder classify"
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
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: (input.size_mb  + 1300) * threads * 2 * attempt
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
  name: "viral-contigident.py PhaMer classify"
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
  threads: 16
  resources:
    mem_mb=lambda wildcards, attempt, threads: attempt * threads * 2 * 10**3
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
  name : "viral-contigident.py merge classification outputs"
  input:
    vs2out=relpath("viralcontigident/samples/{sample_id}/intermediate/vs2/final-viral-score.tsv"),
    dvfout=relpath("viralcontigident/samples/{sample_id}/intermediate/dvf/final_score.txt"), 
    phamerout=relpath("viralcontigident/samples/{sample_id}/intermediate/phamer/out/phamer_prediction.csv"),
  output:
    relpath("viralcontigident/samples/{sample_id}/output/merged_scores.csv")
  params:
    script="workflow/scripts/viralcontigident/mergeout.py",
    vs2minlen=configdict['vs2minlen'],
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
        --vs2out {input.vs2out} \
        --dvfout {input.dvfout} \
        --phamerout {input.phamerout} \
        --vs2minlen {params.vs2minlen} \
        --dvfminlen {params.dvfminlen} \
        --phamerminlen {params.phamerminlen} \
        --output {params.tmpdir}/tmp.csv &> {log}

    mv {params.tmpdir}/* {output}
    rm -rf {params.tmpdir}
    """
    

rule filter_outputs:
  name: "viral-contigident.py filter viral contigs"
  input:
    fna=relpath("viralcontigident/samples/{sample_id}/tmp/final.contigs.filtered.fa"),
    scrs=relpath("viralcontigident/samples/{sample_id}/output/merged_scores.csv")
  output:
    fna=relpath("viralcontigident/samples/{sample_id}/output/viral.contigs.fa"),
    scrs=relpath("viralcontigident/samples/{sample_id}/output/merged_scores_filtered.csv"),
    hits=relpath("viralcontigident/samples/{sample_id}/output/viralhits_list")
  params:
    script="workflow/scripts/viralcontigident/filtercontig_scores_vs2.py",
    vs2_cutoff=configdict['vs2cutoff_p'], 
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
        --vs2_min_score {params.vs2_cutoff} \
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
  name : "viral-contigident.py combine viral contigs"
  input:
    fna=expand(relpath("viralcontigident/samples/{sample_id}/output/viral.contigs.fa"), sample_id=assembly_ids),
    scrs=expand(relpath("viralcontigident/samples/{sample_id}/output/merged_scores_filtered.csv"), sample_id=assembly_ids)
  output: 
    fna=relpath("viralcontigident/intermediate/scores/combined.viralcontigs.fa"),
    scrs=relpath("viralcontigident/intermediate/scores/combined_viral_scores.csv")
  params:
    script="workflow/scripts/viralcontigident/mergeout_scores.py", 
    names=list(assembly_ids),
    tmpdir=tmpd
  log: os.path.join(logdir, "catcontigs.log")
  conda: "../envs/utility.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}

    echo "{params.names}" > {params.tmpdir}/tmp.names
    echo "{input.scrs}" > {params.tmpdir}/tmp.csv.paths

    python {params.script} \
        --names {params.tmpdir}/tmp.names \
        --csvs {params.tmpdir}/tmp.csv.paths > {params.tmpdir}/tmp.csv 2> {log}
    cat {input.fna} > {params.tmpdir}/tmp.fa 2> {log}
    
    mv {params.tmpdir}/tmp.csv {output.scrs}
    mv {params.tmpdir}/tmp.fa {output.fna}

    rm -rf {params.tmpdir}
    """


# THEN GOES INTO 
# 1) VIRAL BINNING [optional]
# 2) CLUSTERING [sensitive or fast]
# 3) CHECKV-PYHMMER 
# THEN COMES BACK HERE TO GET vCONTIGS


rule combine_classifications:
  name: "viral-contigident.py combine derepped classification results"
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
  name: "viral-contigident.py consensus vOTU filtering"
  input:
    relpath("viralcontigident/output/checkv/combined_classification_results.csv")
  output:
    summary=relpath("viralcontigident/output/classification_summary_vOTUs.csv"),
    proviruslist=relpath("viralcontigident/output/provirus.list.txt"),
    viruslist=relpath("viralcontigident/output/virus.list.txt")
  params:
    script="workflow/scripts/viralcontigident/consensus_filtering_vs2.py",
    vs2=configdict['vs2cutoff_s'],
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
        --vs2_min_score {params.vs2} \
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
  name: "viral-contigident.py generate final vOTUs"
  input:
    provirusfasta=relpath("viralcontigident/output/checkv/proviruses.fna"),
    virusfasta=relpath("viralcontigident/output/checkv/viruses.fna"), 
    provirushits=relpath("viralcontigident/output/provirus.list.txt"),
    virushits=relpath("viralcontigident/output/virus.list.txt")
  output:
    combined=relpath("viralcontigident/output/combined.final.vOTUs.fa"),
    provirus=relpath("viralcontigident/output/provirus.final.vOTUs.fa"),
    virus=relpath("viralcontigident/output/virus.final.vOTUs.fa")
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

    rm -rf {params.tmpdir}/*
    """

