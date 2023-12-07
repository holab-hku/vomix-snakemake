configfile: "config/taxonomy.yml"

rule prodigalgv_taxonomy:
  name: "taxonomy.py prodigal-gv vTOUs [parallelized]"
  input: 
    "results/viralcontigident/output/combined.final.vOTUs.fa"
  output: 
    "results/taxonomy/viral/intermediate/prodigal/proteins.vOTUs.faa"
  params:
    script_path="workflow/software/prodigal-gv/parallel-prodigal-gv.py", 
    output_dir="results/taxonomy/viral/intermediate/prodigal/", 
    tmp_dir="$TMPDIR/prodigal"
  conda: "../envs/prodigal-gv.yml"
  log: "logs/taxonomy_prodigalgv.log"
  threads: 64
  resources: 
    mem_mb=lambda wildcards, attempt: attempt * 72 * 10**3
  shell: 
    """
    rm -rf {params.tmp_dir} {params.output_dir}
    mkdir -p {params.tmp_dir} {params.output_dir}

    python {params.script_path} \
        -i {input} \
        -a {params.tmp_dir}/tmp.faa \
        -t {threads} &> {log}

    mv {params.tmp_dir}/tmp.faa {output}

    rm -rf {params.tmp_dir}
    """



rule hmmscan_taxonomy:
  name: "taxonomy.py VIRify hmmscan vOTUs [ViPhOGs]"
  input:
    faa="results/taxonomy/viral/intermediate/prodigal/proteins.vOTUs.faa",
    db="workflow/database/viphogshmm/vpHMM_database_v3/vpHMM_database_v3.hmm"
  output:
    tbl_cutga="results/taxonomy/viral/intermediate/VIRify/vOTUs_hmmscan_cutga.tbl", 
    tbl="results/taxonomy/viral/intermediate/VIRify/vOTUs_hmmscan.tbl"
  params:
    tmp_dir="$TMPDIR/viphogs"
  conda: "../envs/hmmer.yml"
  log: "logs/taxonomy_hmmscan.log"
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 72 * 10**3
  shell:
    """
    rm -rf {params.tmp_dir}
    mkdir -p {params.tmp_dir}
    
    hmmscan --cpu {threads} \
        --noali \
        --cut_ga \
        --domtblout {params.tmp_dir}/tmp_cutga.tbl \
        {input.db} {input.faa} &> {log}

    #filter low evalue>0.001 hits if cut_ga model not available 
    awk "{{if(\$1 ~ /^#/){{print \$0}}else{{if(\$7<0.001){{print \$0}}}}}}" {params.tmp_dir}/tmp_cutga.tbl > tmp.tbl 2> {log}
    
    mv {params.tmp_dir}/tmp_cutga.tbl {output.tbl_cutga}
    mv {params.tmp_dir}/tmp.tbl {output.tbl}

    rm -rf {params.tmp_dir}
    """


rule VIRify_postprocess:
  name: "taxonomy.py VIRify post-process hmmer"
  input: 
    "results/taxonomy/viral/intermediate/VIRify/vOTUs_hmmscan.tbl"
  output:
    "results/taxonomy/viral/intermediate/VIRify/vOTUs_hmmscan_processed.tbl"
  params:
    script_path="workflow/scripts/taxonomy/hmmscan_format_table.py", 
    tmp_dir="$TMPDIR/viphogs"
  conda: "../envs/taxonomy.yml"
  log: "logs/taxonomy_VIRifypostprocess.log"
  shell:
    """
    rm -rf {params.tmp_dir}
    mkdir -p {params.tmp_dir}

    python {params.script_path} -t {input} -o {params.tmp_dir}/tmp.tbl &> {log}
    mv {params.tmp_dir}/tmp.tbl {output}

    rm -rf {params.tmp_dir}
    """

rule VIRify_ratioeval:
  name: "taxonomy.py VIRify calculate evalue ratio"
  input:
    tbl="results/taxonomy/viral/intermediate/VIRify/vOTUs_hmmscan_processed.tbl", 
    tsv="workflow/database/viphogshmm/additional_data_vpHMMs_v4.tsv" 
  output:
    "results/taxonomy/viral/intermediate/VIRify/vOTUs_hmmscan_processed_modified_informative.tsv"
  params:
    script_path="workflow/scripts/taxonomy/ratio_evalue_table.py", 
    tmp_dir="$TMPDIR/viphogs",
    evalue=config['VIRifyhmmeval']
  conda: "../envs/taxonomy.yml"
  log: "logs/taxonomy_VIRifyratioeval.log"
  shell:
    """
    rm -rf {params.tmp_dir}
    mkdir -p {params.tmp_dir}

    python {params.script_path} \
        -i {input.tbl} \
        -t {input.tsv} \
        -e {params.evalue} \
        -o {params.tmp_dir}/tmp.tsv &> {log}
    mv {params.tmp_dir}/tmp.tsv {output}

    rm -rf {params.tmp_dir}
    """


rule VIRify_annotate:
  name: "taxonomy.py VIRify annotate proteins"
  input:
    faa="results/taxonomy/viral/intermediate/prodigal/proteins.vOTUs.faa", 
    tsv="results/taxonomy/viral/intermediate/VIRify/vOTUs_hmmscan_processed_modified_informative.tsv"
  output:
    "results/taxonomy/viral/intermediate/VIRify/vOTUs_annotation.tsv"
  params:
    script_path="workflow/scripts/taxonomy/viral_contigs_annotation.py",
    tmp_dir="$TMPDIR/viphogs"
  conda: "../envs/taxonomy.yml"
  log: "logs/taxonomy_VIRifyannotation.log"
  shell:
    """
    rm -rf {params.tmp_dir}
    mkdir -p {params.tmp_dir}

    python {params.script_path} \
        -p {input.faa} \
        -t {input.tsv} \
        -o {params.tmp_dir}/tmp.tsv
    mv {params.tmp_dir}/tmp.tsv {output}

    rm -rf {params.tmp_dir}
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
    script_path="workflow/scripts/taxonomy/contig_taxonomic_assign.py",
    output_dir="results/taxonomy/viral/output/VIRify/",
    thresh=config['VIRifyprop'], 
    tmp_dir="$TMPDIR/viphogs"
  conda: "../envs/taxonomy.yml"
  log: "logs/taxonomy_VIRifyassign.log"
  shell:
    """
    rm -r {params.tmp_dir} {params.output_dir}
    mkdir -p {params.tmp_dir} {params.output_dir}

    python {params.script_path} \
        -i {input.tsv} \
        -d {input.db} \
        --factor {input.csv} \
        --taxthres {params.thresh} \
        -o {params.tmp_dir}/tmp.tsv 2> {log}
    mv {params.tmp_dir}/tmp.tsv {output}

    rm -r {params.tmp_dir}
    """

# rule blastp_taxonomy:
#   input: 
#     "results/viralcontigident/output/combined.final.vOTUs.fa"
#   output: 
#     "results/viralcontigident/output/combined.final.vOTUs.fa"
#   params:
#     tmp_dir="$TMPDIR"
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
