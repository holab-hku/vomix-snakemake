import os 
configfile: "config/viralcontigident.yml"

logdir = relpath("viralcontigident/logs")
benchmarks = relpath("viralcontigident/benchmarks")
tmpd = relpath("viralcontigident/tmp")

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

############################
# Single-Sample Processing #
############################

if config["contigfile"] != "":
  path = config["contigfile"]
  cwd = os.getcwd()
  path_full = os.path.join(cwd, path)

  if not os.path.exists(path):
    print("The contig file path provided does not exist.")
    sys.exit(1)

  sample_id = os.path.basename(path).rsplit(".", 1)[0]
  outdir_p = os.path.join(cwd, relpath("assembly/"), sample_id, "output/")
  outfile = os.path.join(outdir_p, "final.contigs.fa")

  os.makedirs(outdir_p, exist_ok=True)

  try:
    os.symlink(path_full, outfile)
  except FileExistsError:
    print(f"Contig path already exists '{outdir}final.contigs.fa'. The existing file will be used for viral contig annotation. If needed, please delete the current contig file path and try again.")
    pass

  sample_ids = [sample_id]

else:
  sample_ids = samples.keys()
  
###########################
# Multi-sample Processing #
###########################

rule filter_contigs:
  name: "viralcontigident.py filter contigs [length]"
  localrule: True
  input:
    relpath("assembly/{sample_id}/output/final.contigs.fa")
  output:
    relpath("viralcontigident/samples/{sample_id}/tmp/final.contigs.filtered.fa")
  params:
    minlen=config['contigminlen'],
    outdir=relpath("viralcontigident/samples/{sample_id}/tmp"), 
    tmpdir=os.path.join(tmpd, "{sample_id}")
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
  name: "viralcontigident.py geNomad classify" 
  input:
    fna=relpath("assembly/{sample_id}/output/final.contigs.fa")
    #fna=relpath("viralcontigident/samples/{sample_id}/tmp/final.contigs.filtered.fa"),
  output:
    relpath("viralcontigident/samples/{sample_id}/intermediate/genomad/final.contigs_summary/final.contigs_virus_summary.tsv")
  params:
    genomadparams=config['genomadparams'],
    dbdir=config['genomaddb'],
    outdir=relpath("viralcontigident/samples/{sample_id}/intermediate/genomad/"),
    tmpdir=os.path.join(tmpd, "{sample_id}")
  log: os.path.join(logdir, "genomad_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "genomad_{sample_id}.log")
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
    fna=relpath("assembly/{sample_id}/output/final.contigs.fa"), 
    #fna=relpath("viralcontigident/samples/{sample_id}/tmp/final.contigs.filtered.fa")
  output:
    relpath("viralcontigident/samples/{sample_id}/intermediate/dvf/final_score.txt")
  params:
    script="workflow/software/DeepVirFinder/dvf.py",
    minlen=config['contigminlen'],
    parameters=config['dvfparams'], 
    modeldir="workflow/software/DeepVirFinder/models/",
    outdir=relpath("viralcontigident/samples/{sample_id}/intermediate/dvf/"),
    tmpdir=os.path.join(tmpd, "{sample_id}")
  log: os.path.join(logdir, "dvf_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "dvf_{sample_id}.log")
  conda: "../envs/dvf.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
  shell:
    """
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        -i {input.fna} \
        -l {params.minlen} \
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
    fna=relpath("assembly/{sample_id}/output/final.contigs.fa")
    #fna=relpath("viralcontigident/samples/{sample_id}/tmp/final.contigs.filtered.fa")
  output:
    relpath("viralcontigident/samples/{sample_id}/intermediate/phamer/out/phamer_prediction.csv")
  params:
    script="workflow/software/PhaBOX/PhaMer_single.py",
    scriptdir="workflow/software/PhaBOX/scripts/",
    parameters=config['phamerparams'],
    paramsdir="workflow/params/phabox/",
    minlen=config['contigminlen'],
    dbdir=config['phamerdb'],
    outdir=relpath("viralcontigident/samples/{sample_id}/intermediate/phamer/"),
    tmpdir=os.path.join(tmpd, "{sample_id}")
  log: os.path.join(logdir, "phamer_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "phamer_{sample_id}.log")
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
        --len {params.minlen} \
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
    genomadout=relpath("viralcontigident/samples/{sample_id}/intermediate/genomad/final.contigs_summary/final.contigs_virus_summary.tsv"),
    dvfout=relpath("viralcontigident/samples/{sample_id}/intermediate/dvf/final_score.txt"), 
    phamerout=relpath("viralcontigident/samples/{sample_id}/intermediate/phamer/out/phamer_prediction.csv"),
  output:
    relpath("viralcontigident/samples/{sample_id}/output/merged_scores.csv")
  params:
    script="workflow/scripts/viralcontigident/mergeout.py",
    genomadminlen=config['genomadminlen'],
    dvfminlen=config['dvfminlen'], 
    phamerminlen=config['dvfminlen'], 
    outdir=relpath("viralcontigident/samples/{sample_id}/output/"),
    tmpdir=os.path.join(tmpd, "{sample_id}")
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
  name : "viralcontigident.py filter viral contigs"
  input:
    fna=relpath("viralcontigident/samples/{sample_id}/tmp/final.contigs.filtered.fa"),
    scrs=relpath("viralcontigident/samples/{sample_id}/output/merged_scores.csv")
  output:
    fna=relpath("viralcontigident/samples/{sample_id}/output/viral.contigs.fa"),
    scrs=relpath("viralcontigident/samples/{sample_id}/output/merged_scores_filtered.csv"),
    hits=relpath("viralcontigident/samples/{sample_id}/output/viralhits_list")
  params:
    script="workflow/scripts/viralcontigident/filtercontig_scores.py",
    genomad_cutoff=config['genomadcutoff'], 
    dvf_cutoff=config['dvfcutoff'], 
    dvf_pvalmax=config['dvfpval'],
    phamer_cutoff=config['phamercutoff'], 
    phamer_pred=config['phamerpred'], 
    tmpdir=os.path.join(tmpd, "{sample_id}")
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
  name : "viralcontigident.py combine viral contigs"
  input:
    fna=expand(relpath("viralcontigident/samples/{sample_id}/output/viral.contigs.fa"), sample_id = sample_ids),
    scrs=expand(relpath("viralcontigident/samples/{sample_id}/output/merged_scores_filtered.csv"), sample_id = sample_ids)
  output: 
    fna=relpath("viralcontigident/intermediate/scores/combined.viralcontigs.fa"),
    scrs=relpath("viralcontigident/intermediate/scores/combined_viral_scores.csv")
  params:
    script="workflow/scripts/viralcontigident/mergeout_scores.py", 
    tmpdir=tmpd
  log: os.path.join(logdir, "catcontigs.log")
  conda: "../envs/utility.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}

    python {params.script} {input.scrs} > {params.tmpdir}/tmp.csv 2> {log}
    cat {input.fna} > {params.tmpdir}/tmp.fa 2> {log}

    mv {params.tmpdir}/tmp.csv {output.scrs}
    mv {params.tmpdir}/tmp.fa {output.fna}

    rm -rf {params.tmpdir}
    """


rule checkv:
  name : "viralcontigident.py CheckV dereplicated contigs"
  input:
    relpath("viralcontigident/output/derep/combined.viralcontigs.derep.fa")
  output:
    relpath("viralcontigident/output/checkv/viruses.fna"),
    relpath("viralcontigident/output/checkv/proviruses.fna"), 
    relpath("viralcontigident/output/checkv/quality_summary.tsv")
  params:
    checkvparams= config['checkvparams'],
    outdir=relpath("viralcontigident/output/checkv"),
    tmpdir=os.path.join(tmpd, "checkv"),
    dbdir="workflow/database/checkv"
  log: os.path.join(logdir, "checkv.log")
  benchmark: os.path.join(benchmarks, "checkv.log")
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
  name: "viralcontigident.py consensus vOTU filtering"
  input:
    relpath("viralcontigident/output/checkv/combined_classification_results.csv")
  output:
    summary=relpath("viralcontigident/output/classification_summary_vOTUs.csv"),
    proviruslist=relpath("viralcontigident/output/provirus.list.txt"),
    viruslist=relpath("viralcontigident/output/virus.list.txt")
  params:
    script="workflow/scripts/viralcontigident/consensus_filtering.py",
    genomad=config['genomadcutoff'],
    dvf=config['dvfcutoff'],
    phamer=config['phamerpred'], 
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
        --phamer_pred {params.phamer} \
        --summary_out {params.tmpdir}/tmp.csv \
        --provirus_list {params.tmpdir}/tmp.list.1 \
        --virus_list {params.tmpdir}/tmp.list.2  2> {log}

    mv {params.tmpdir}/tmp.csv {output.summary}
    mv {params.tmpdir}/tmp.list.1 {output.proviruslist}
    mv {params.tmpdir}/tmp.list.2 {output.viruslist}
    """


rule votu:
  name: "viralcontigident.py generate final vOTUs"
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

    seqkit replace {input.provirusfasta} --f-use-regexp -p "(.+)_\d\s.+$" -r '$1' | \
        seqkit grep -f {input.provirushits} > {params.tmpdir}/tmp1.fa 2> {log}
    seqkit replace {input.virusfasta} --f-use-regexp -p "(.+)_\d\s.+$" -r '$1' | \
        seqkit grep -f {input.provirushits} >> {params.tmpdir}/tmp1.fa 2> {log}

    seqkit grep {input.virusfasta} -f {input.virushits} > {params.tmpdir}/tmp2.fa 2> {log}
    seqkit grep {input.provirusfasta} -f {input.virushits} >> {params.tmpdir}/tmp2.fa 2> {log}

    cat {params.tmpdir}/tmp1.fa {params.tmpdir}/tmp2.fa  > {params.tmpdir}/tmp3.fa 2> {log}

    seqkit rmdup {params.tmpdir}/tmp1.fa > {output.provirus} 2> {log}
    seqkit rmdup {params.tmpdir}/tmp2.fa > {output.virus} 2> {log}
    seqkit rmdup {params.tmpdir}/tmp3.fa > {output.combined} 2> {log}

    rm -rf {params.tmpdir}/*
    """

rule done_log:
  name: "viralcontigident.py removing tmp files"
  input:
    pseudo=relpath("viralcontigident/output/combined.final.vOTUs.fa")
  output:
    os.path.join(logdir, "done.log")
  params:
    filteredcontigs=expand(relpath("viralcontigident/samples/{sample_id}/tmp"), sample_id = sample_ids),
    tmpdir=tmpd
  log: os.path.join(logdir, "done.log")
  shell:
    """
    rm -rf {params.filteredcontigs} {params.tmpdir} /*
    touch {output}
    """

