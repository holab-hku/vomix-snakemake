configfile: "config/taxonomy.yml"
logdir = relpath("taxonomy/viral/logs")
tmpd = relpath("taxonomy/viral/tmp")
os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

rule prodigalgv_taxonomy:
  name: "taxonomy.py prodigal-gv vTOUs [parallelized]"
  input: 
    relpath("viralcontigident/output/combined.final.vOTUs.fa")
  output: 
    relpath("taxonomy/viral/intermediate/prodigal/proteins.vOTUs.faa")
  params:
    script="workflow/software/prodigal-gv/parallel-prodigal-gv.py", 
    outdir=relpath("taxonomy/viral/intermediate/prodigal/"),
    tmpdir=os.path.join(tmpd, "prodigal")
  conda: "../envs/prodigal-gv.yml"
  log: os.path.join(logdir, "prodigal-gv.log")
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
    faa=relpath("taxonomy/viral/intermediate/prodigal/proteins.vOTUs.faa"),
    db="workflow/database/viphogshmm/vpHMM_database_v3/vpHMM_database_v3.hmm"
  output:
    tblcutga=relpath("taxonomy/viral/intermediate/viphogs/vOTUs_hmmsearch_cutga.tbl"), 
    tbl=relpath("taxonomy/viral/intermediate/viphogs/vOTUs_hmmscan.tbl")
  params:
    script="workflow/scripts/taxonomy/pyhmmer_wrapper.py", 
    outdir=relpath("taxonomy/viral/intermediate/viphogs/"),
    tmpdir=os.path.join(tmpd, "viphogs")
  conda: "../envs/pyhmmer.yml"
  log: os.path.join(logdir, "pyhmmer.log")
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
    relpath("taxonomy/viral/intermediate/viphogs/vOTUs_hmmscan.tbl")
  output:
    relpath("taxonomy/viral/intermediate/viphogs/vOTUs_hmmsearch_processed.tsv")
  params:
    script="workflow/scripts/taxonomy/hmmscan_format_table.py", 
    tmpdir=os.path.join(tmpd, "viphogs")
  conda: "../envs/taxonomy.yml"
  log: os.path.join(logdir, "VIRify_postprocess.log")
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
    tbl=relpath("taxonomy/viral/intermediate/viphogs/vOTUs_hmmsearch_processed.tsv"), 
    tsv="workflow/database/viphogshmm/additional_data_vpHMMs_v4.tsv" 
  output:
    relpath("taxonomy/viral/intermediate/viphogs/vOTUs_hmmsearch_processed_modified_informative.tsv")
  params:
    script="workflow/scripts/taxonomy/ratio_evalue_table.py", 
    tmpdir=os.path.join(tmpd, "viphogs"),
    evalue=config['viphogshmmeval']
  conda: "../envs/taxonomy.yml"
  log: os.path.join(logdir, "VIRify_ratioeval.log")
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
    faa=relpath("taxonomy/viral/intermediate/prodigal/proteins.vOTUs.faa"),
    tsv=relpath("taxonomy/viral/intermediate/viphogs/vOTUs_hmmsearch_processed_modified_informative.tsv")
  output:
    relpath("taxonomy/viral/intermediate/viphogs/vOTUs_annotation.tsv")
  params:
    script="workflow/scripts/taxonomy/viral_contigs_annotation.py",
    tmpdir=os.path.join(tmpd, "viphogs")
  conda: "../envs/taxonomy.yml"
  log: os.path.join(logdir, "VIRify_annotation.log")
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
    tsv=relpath("taxonomy/viral/intermediate/viphogs/vOTUs_annotation.tsv"), 
    db="workflow/database/ncbi/ete3ncbitax.sqlite", 
    csv="workflow/params/VIRify/viphogs_cds_per_taxon_cummulative.csv"
  output:
    relpath("taxonomy/viral/intermediate/viphogs/taxonomy.tsv")
  params:
    script="workflow/scripts/taxonomy/contig_taxonomic_assign.py",
    outdir=relpath("taxonomy/viral/intermediate/viphogs/"),
    thresh=config['viphogsprop'], 
    tmpdir=os.path.join(tmpd, "viphogs")
  conda: "../envs/taxonomy.yml"
  log: os.path.join(logdir, "VIRify_assign.log")
  shell:
    """
    rm -rf {params.tmpdir} {output}
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


rule genomad_taxonomy:
  name: "taxonomy.py Parse geNomad taxonomy"
  input:
    relpath("viralcontigident/output/classification_summary_vOTUs.csv")
  output:
    relpath("taxonomy/viral/intermediate/genomad/taxonomy.csv")
  params:
    script="workflow/scripts/taxonomy/genomad_taxonomy.py",
    outdir=relpath("taxonomy/viral/intermediate/genomad"),
    tmpdir=os.path.join(tmpd, "genomad")
  conda: "../envs/taxonomy.yml"
  log: os.path.join(logdir, "genomad_parse.log")
  shell:
    """
    rm -rf {params.tmpdir}/* {params.outdir}/*
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        --input {input} \
        --output {params.tmpdir}/tmp.csv 2> {log}

    mv {params.tmpdir}/tmp.csv {output}
    rm -rf {params.tmpdir}/*

    """


rule phagcn_taxonomy:
  name: "taxonomy.py PhaGCN phage taxonomy"
  input:
    fna=relpath("viralcontigident/output/combined.final.vOTUs.fa")
  output:
    relpath("taxonomy/viral/intermediate/phagcn/taxonomy.csv")
  params:
    script="workflow/software/PhaBOX/PhaGCN_single.py",
    scriptdir="workflow/software/PhaBOX/scripts/",
    parameters=config['phagcnparams'],
    paramsdir="workflow/params/phabox/",
    dbdir=config['phagcndb'],
    outdir=relpath("taxonomy/viral/intermediate/phagcn"),
    tmpdir=os.path.join(tmpd, "phagcn")
  log: os.path.join(logdir, "phagcn.log")
  conda: "../envs/phabox.yml"
  threads: 16
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 32 * 10**3
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        --contigs {input.fna} \
        --threads {threads} \
        --rootpth {params.tmpdir} \
        --dbdir {params.dbdir} \
        --parampth {params.paramsdir} \
        --scriptpth {params.scriptdir} \
        {params.parameters} &> {log}

    mv {params.tmpdir}/* {params.outdir}/
    cp {params.outdir}/out/phagcn_prediction.csv {output}

    rm -rf {params.tmpdir}
    """ 



rule diamond_makedb:
  name: "taxonomy.py make NCBI-Virus Refseq proteins diamond database"
  input:
    "workflow/database/ncbi/ncbi-virus/ncbi.virus.RefSeq.faa"
  output:
    "workflow/database/diamond/ncbi.virus.Refseq.dmnd"
  params:
    outdir="workflow/database/diamond",
    tmpdir=os.path.join(tmpd, "diamond")
  log: os.path.join(logdir, "diamond_makedb.log")
  conda: "../envs/diamond.yml"
  threads: 32
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    diamond makedb --in {input} --db {params.tmpdir}/tmp.dmnd --threads {threads}

    mv {params.tmpdir}/tmp.dmnd {output}
    """

rule dimaond_taxonomy:
  name: "taxonomy.py DIAMOND blastp [NCBI-Virus Refseq]"
  input:
    faa=relpath("taxonomy/viral/intermediate/prodigal/proteins.vOTUs.faa"),
    db="workflow/database/diamond/ncbi.virus.Refseq.dmnd"
  output:
    relpath("taxonomy/viral/intermediate/diamond/diamond_out.tsv")
  params:
    parameters=config['diamondparams'],
    outdir=relpath("taxonomy/viral/intermediate/diamond"),
    tmpdir=os.path.join(tmpd, "diamond")
  log: os.path.join(logdir, "diamond_blastp.log")
  conda: "../envs/diamond.yml"
  threads: 32
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 32 * 10**3
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    diamond blastp \
        --header \
        --db {input.db} \
        --query {input.faa} \
        --out {params.tmpdir}/tmp.tsv \
        --threads {threads} \
        {params.parameters} &> {log}

    mv {params.tmpdir}/tmp.tsv {output}
    """



rule diamond_parse_taxonomy:
  name: "taxonomy.py DIAMOND taxonomic classification [NCBI-Virus Refseq]"
  input:
    diamondout = relpath("taxonomy/viral/intermediate/diamond/diamond_out.tsv"),
    taxcsv = "workflow/database/ncbi/ncbi-virus/refseq/refseq.metadata.csv"
  output:
    relpath("taxonomy/viral/intermediate/diamond/taxonomy.csv")
  params:
    script="workflow/scripts/taxonomy/diamond_taxonomy_parse.py",
    outdir=relpath("taxonomy/viral/intermediate/diamond"),
    thresh=0.7, 
    taxcols="species,genus,family,order,class,phylum,kingdom",
    tmpdir=os.path.join(tmpd, "diamond")
  log: os.path.join(logdir, "diamond_assign_taxonomy.log")
  conda: "../envs/taxonomy.yml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 32 * 10**3
  shell:
    """
    rm -rf {params.tmpdir} {output}
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        --diamondout {input.diamondout} \
        --taxcsv {input.taxcsv} \
        --threshold {params.thresh} \
        --taxcolumns {params.taxcols} \
        --outputcsv {params.tmpdir}/tmp.csv &> {log}

    mv {params.tmpdir}/tmp.csv {output}
    """



rule merge_taxonomy:
  name: "taxonomy.py merge taxonomic classifications"
  localrule: True
  input:
    diamond=relpath("taxonomy/viral/intermediate/diamond/taxonomy.csv"),
    viphogs=relpath("taxonomy/viral/intermediate/viphogs/taxonomy.tsv"),
    phagcn=relpath("taxonomy/viral/intermediate/phagcn/taxonomy.csv"),
    genomad=relpath("taxonomy/viral/intermediate/genomad/taxonomy.csv"),
    contigs=relpath("viralcontigident/output/combined.final.vOTUs.fa")
  output:
    relpath("taxonomy/viral/output/merged_taxonomy.csv")
  params:
    script="workflow/scripts/taxonomy/merge_taxonomy.py",
    outdir=relpath("taxonomy/viral/output/"),
    tmpdir=os.path.join(tmpd, "merge")
  log: os.path.join(logdir, "merge_taxonomy.log")
  conda: "../envs/taxonomy.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir}/* {output}
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        --diamondout {input.diamond} \
        --viphogsout {input.viphogs} \
        --phagcnout {input.phagcn} \
        --genomadout {input.genomad} \
        --contigs {input.contigs} \
        --output {params.tmpdir}/tmp.csv
    
    mv {params.tmpdir}/tmp.csv {output}
    """
        


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
