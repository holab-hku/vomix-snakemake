configfile: "config/taxonomy.yml"

rule prodigalgv_taxonomy:
  name: "taxonomy.py prodigal-gv vTOUs [parallelized]"
  input: 
    "results/viralcontigident/output/combined.final.vOTUs.fa"
  output: 
    "results/taxonomy/viral/intermediate/prodigal/proteins.vOTUs.faa"
  params:
    script="workflow/software/prodigal-gv/parallel-prodigal-gv.py", 
    outdir="results/taxonomy/viral/intermediate/prodigal/", 
    tmpdir="$TMPDIR/prodigal"
  conda: "../envs/prodigal-gv.yml"
  log: "logs/taxonomy_prodigalgv.log"
  threads: 64
  resources: 
    mem_mb=lambda wildcards, attempt: attempt * 72 * 10**3
  shell: 
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        -i {input} \
        -a {params.tmpdir}/tmp.faa \
        -t {threads} &> {log}

    mv {params.tmpdir}/tmp.faa {output}

    rm -rf {params.tmpdir}
    """



rule pyhmmer_taxonomy:
  name: "taxonomy.py PyHMMER vOTUs [ViPhOGs]"
  input:
    faa="results/taxonomy/viral/intermediate/prodigal/proteins.vOTUs.faa",
    db="workflow/database/viphogshmm/vpHMM_database_v3/vpHMM_database_v3.hmm"
  output:
    tblcutga="results/taxonomy/viral/intermediate/VIRify/vOTUs_hmmsearch_cutga.tbl", 
    tbl="results/taxonomy/viral/intermediate/VIRify/vOTUs_hmmscan.tbl"
  params:
    script="workflow/scripts/taxonomy/pyhmmer_wrapper.py", 
    outdir="results/taxonomy/viral/intermediate/VIRify/", 
    tmpdir="$TMPDIR/viphogs"
  conda: "../envs/pyhmmer.yml"
  log: "logs/taxonomy_pyhmmer.log"
  threads: 32
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 72 * 10**3
  shell:
    """
    rm -rf {params.tmpdir} {output.tblcutga} {output.tbl}
    mkdir -p {params.tmpdir}
    
    python {params.script} \
        --proteins {input.faa} \
        --hmmdb {input.db} \
        --bit_cutoff gathering \
        --cores {threads} \
        --hmmscan \
        --domtblout {params.tmpdir}/tmpcutga.tbl &> {log}

    #filter low evalue>0.001 hits if cut_ga model not available 
    awk "{{if(\$1 ~ /^#/){{print \$0}}else{{if(\$7<0.001){{print \$0}}}}}}" {params.tmpdir}/tmpcutga.tbl > {params.tmpdir}/tmp.tbl 2>> {log}
    
    mv {params.tmpdir}/tmpcutga.tbl {output.tblcutga}
    mv {params.tmpdir}/tmp.tbl {output.tbl}

    rm -rf {params.tmpdir}
    """


rule VIRify_postprocess:
  name: "taxonomy.py VIRify post-process hmmer"
  input: 
    "results/taxonomy/viral/intermediate/VIRify/vOTUs_hmmscan.tbl"
  output:
    "results/taxonomy/viral/intermediate/VIRify/vOTUs_hmmsearch_processed.tsv"
  params:
    script="workflow/scripts/taxonomy/hmmscan_format_table.py", 
    tmpdir="$TMPDIR/viphogs"
  conda: "../envs/taxonomy.yml"
  log: "logs/taxonomy_VIRifypostprocess.log"
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}

    python {params.script} -t {input} -o {params.tmpdir}/tmp.tbl &> {log}
    mv {params.tmpdir}/tmp.tbl {output}

    rm -rf {params.tmpdir}
    """

rule VIRify_ratioeval:
  name: "taxonomy.py VIRify calculate evalue ratio"
  input:
    tbl="results/taxonomy/viral/intermediate/VIRify/vOTUs_hmmsearch_processed.tsv", 
    tsv="workflow/database/viphogshmm/additional_data_vpHMMs_v4.tsv" 
  output:
    "results/taxonomy/viral/intermediate/VIRify/vOTUs_hmmsearch_processed_modified_informative.tsv"
  params:
    script="workflow/scripts/taxonomy/ratio_evalue_table.py", 
    tmpdir="$TMPDIR/viphogs",
    evalue=config['VIRifyhmmeval']
  conda: "../envs/taxonomy.yml"
  log: "logs/taxonomy_VIRifyratioeval.log"
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}

    python {params.script} \
        -i {input.tbl} \
        -t {input.tsv} \
        -e {params.evalue} \
        -o {params.tmpdir}/tmp.tsv &> {log}
    mv {params.tmpdir}/tmp.tsv {output}

    rm -rf {params.tmpdir}
    """


rule VIRify_annotate:
  name: "taxonomy.py VIRify annotate proteins"
  input:
    faa="results/taxonomy/viral/intermediate/prodigal/proteins.vOTUs.faa", 
    tsv="results/taxonomy/viral/intermediate/VIRify/vOTUs_hmmsearch_processed_modified_informative.tsv"
  output:
    "results/taxonomy/viral/intermediate/VIRify/vOTUs_annotation.tsv"
  params:
    script="workflow/scripts/taxonomy/viral_contigs_annotation.py",
    tmpdir="$TMPDIR/viphogs"
  conda: "../envs/taxonomy.yml"
  log: "logs/taxonomy_VIRifyannotation.log"
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}

    python {params.script} \
        -p {input.faa} \
        -t {input.tsv} \
        -o {params.tmpdir}/tmp.tsv
    mv {params.tmpdir}/tmp.tsv {output}

    rm -rf {params.tmpdir}
    """
    
rule VIRify_assign:
  name: "taxonomy.py VIRify assign taxonomy"
  input:
    tsv="results/taxonomy/viral/intermediate/VIRify/vOTUs_annotation.tsv", 
    db="workflow/database/ncbi/ete3ncbitax.sqlite", 
    csv="workflow/params/VIRify/viphogs_cds_per_taxon_cummulative.csv"
  output:
    "results/taxonomy/viral/output/VIRify/taxonomy.tsv"
  params:
    script="workflow/scripts/taxonomy/contig_taxonomic_assign.py",
    outdir="results/taxonomy/viral/output/VIRify/",
    thresh=config['VIRifyprop'], 
    tmpdir="$TMPDIR/viphogs"
  conda: "../envs/taxonomy.yml"
  log: "logs/taxonomy_VIRifyassign.log"
  shell:
    """
    rm -r {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        -i {input.tsv} \
        -d {input.db} \
        --factor {input.csv} \
        --taxthres {params.thresh} \
        -o {params.tmpdir}/tmp.tsv 2> {log}
    mv {params.tmpdir}/tmp.tsv {output}

    rm -r {params.tmpdir}
    """

# rule blastp_taxonomy:
#   input: 
#     "results/viralcontigident/output/combined.final.vOTUs.fa"
#   output: 
#     "results/viralcontigident/output/combined.final.vOTUs.fa"
#   params:
#     tmpdir="$TMPDIR"
#   conda: "../envs/blast.yml"
#   log:
#   threads:
#   shell:
#     """
#     """
# 
# 
# rule vcontact2_taxonomy:
#   input:
#   output:
#   params:
#   conda: "../envs/vcontact2.yml"
#   log:
#   threads:
#   shell:
#     """
#     """
# 
# rule viphogs_taxonomy:
#   input:
#   output:
#   params:
#   conda: "../envs/viphogs.yml"
#   log:
#   threads:
#   shell:
#     """
#     """
# 
# rule phagcn_taxonomy:
#   input:
#   output:
#   params:
#   conda: "../envs/phabox.yml"
#   log:
#   threads:
#   shell:
#     """
#     """
# 
# rule merge_taxonomy:
#   input:
#   output:
#   params:
#   conda:
#   log:
#   threads:
#   shell:
#     """
#     """
# 
# rule consensus_taxonomy:
#   input:
#   output:
#   params:
#   conda:
#   log:
#   threads:
#   shell:
#     """
#     """
