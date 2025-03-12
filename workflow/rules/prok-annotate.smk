logdir=relpath("annotate/prok/logs")
benchmarks=relpath("annotate/prok/benchmarks")
tmpd = relpath("annotate/prok/tmp")

email=config["email"]
nowstr=config["latest_run"]
outdir=config["outdir"] 
datadir=config["datadir"]

samples, assemblies = parse_sample_list(config["samplelist"], datadir, outdir, email, nowstr)

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)


### MASTER RULE 
rule done_log:
  name: "prok-annotate.smk Done. removing tmp files"
  localrule: True
  input:
    expand(relpath("annotate/prok/samples/{sample_id}/input/{sample_id}_merged.fastq.gz"), sample_id = samples.keys()), 
    expand(relpath("annotate/prok/samples/{sample_id}/{sample_id}_{output_type}.tsv"), sample_id = samples.keys(), output_type = ["pathabundance", "pathcoverage", "genefamilies"]), 
    expand(relpath("annotate/prok/output/{output_type}_merged.tsv"), output_type = ["pathabundance", "pathcoverage", "genefamilies"]),
    expand(relpath("annotate/prok/output/genefamilies_merged-cpm-{database}.tsv"), database = ["ec", "eggnog", "go", "ko", "pfam"]),
    expand(relpath("annotate/prok/output/genefamilies_merged-cpm-{database}-named.tsv"), database = ["ec", "eggnog", "go", "ko", "pfam"]),
  output:
    os.path.join(logdir, "done.log")
  params:
    tmpdir=tmpd
  log: os.path.join(logdir, "done.log")
  shell:
    """
    rm -rf {params.tmpdir}/*
    touch {output}
    """


### RULES
rule pair_fastq:
  name: "prok-annotate.smk merge paired files"
  input: 
    R1=relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"),
    R2=relpath("preprocess/samples/{sample_id}/{sample_id}_R2.fastq.gz"),
  output:
    relpath("annotate/prok/samples/{sample_id}/input/{sample_id}_merged.fastq.gz")
  params:
    outdir=relpath("annotate/prok/samples/{sample_id}/input"),
    tmpdir=os.path.join(tmpd, "cat/{sample_id}")
  log: os.path.join(logdir, "paired_concat_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "paired_concat_{sample_id}.log")
  threads: 1
  resources:
    mem_mb = lambda wildcards, input: max(4*input.size_mb, 2000)
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}
    mkdir -p {params.outdir} {params.tmpdir}

    cat {input.R1} {input.R2} > {params.tmpdir}/tmp.fastq.gz 2> {log}
    mv {params.tmpdir}/tmp.fastq.gz {output}
    """


rule humann3:
  name: "prok-annotate.smk HUMAnN3 run."
  input:
    fastq=relpath("annotate/prok/samples/{sample_id}/input/{sample_id}_merged.fastq.gz"),
    chocophlandb=os.path.join(config['humann-db'], "chocophlan/alaS.centroids.v201901_v31.ffn.gz"), 
    unirefdb=os.path.join(config['humann-db'], "uniref/uniref90_201901b_full.dmnd"), 
    utilitydb=os.path.join(config['humann-db'], "utility_mapping/map_ec_name.txt.gz")
  output:
    gene=relpath("annotate/prok/samples/{sample_id}/{sample_id}_genefamilies.tsv"),
    abund=relpath("annotate/prok/samples/{sample_id}/{sample_id}_pathabundance.tsv"),
    cov=relpath("annotate/prok/samples/{sample_id}/{sample_id}_pathcoverage.tsv")
  params:
    parameters=config["humann-params"],
    outdir=relpath("annotate/prok/samples/{sample_id}"), 
    dbdir=config['humann-db'], 
    name="{sample_id}",
    tmpdir=os.path.join(tmpd, "humann/{sample_id}"), 
  conda: "../envs/biobakery3.yml"
  log: os.path.join(logdir, "HUMAnN3_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "HUMAnN3_{sample_id}.log")
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt, threads, input: max(64 * 10**3 + 4000, 4000)
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir} {params.outdir}

    humann \
        --input {input.fastq} \
        --output {params.tmpdir} \
        --threads {threads} \
        --output-format tsv \
        --output-basename {params.name} \
        --protein-database {params.dbdir}/uniref \
        --nucleotide-database {params.dbdir}/chocophlan \
        --log-level INFO \
        --o-log {log} \
        {params.parameters} &>> {log}

    mv {params.tmpdir}/* {params.outdir}/
    rm -rf {params.tmpdir}
    """


rule humann_merge:
  name: "prok-annotate.smk HUMAnN3 merge results"
  input:
    gene=expand(relpath("annotate/prok/samples/{sample_id}/{sample_id}_genefamilies.tsv"), sample_id = assemblies.keys()),
    abund=expand(relpath("annotate/prok/samples/{sample_id}/{sample_id}_pathabundance.tsv"), sample_id = assemblies.keys()),
    cov=expand(relpath("annotate/prok/samples/{sample_id}/{sample_id}_pathcoverage.tsv"), sample_id = assemblies.keys())
  output:
    gene=relpath("annotate/prok/output/genefamilies_merged.tsv"),
    abund=relpath("annotate/prok/output/pathabundance_merged.tsv"),
    cov=relpath("annotate/prok/output/pathcoverage_merged.tsv")
  params:
    outdir=relpath("annotate/prok/output"),
    tmpdir=os.path.join(tmpd, "humann_merge")
  conda: "../envs/biobakery3.yml"
  log: os.path.join(logdir, "humann_merge.log")
  benchmark: os.path.join(benchmarks, "humann_merge.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, threads, input: input.size_mb + 4000
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}
    mkdir -p {params.outdir} {params.tmpdir}
  
    for file in {input.gene}; do ln -s $(pwd)/$file $(pwd)/{params.tmpdir} &>> {log}; done
    humann_join_tables --input {params.tmpdir} --output {output.gene} &>> {log}
    rm {params.tmpdir}/* 

    for file in {input.abund}; do ln -s $(pwd)/$file $(pwd)/{params.tmpdir} &>> {log}; done
    humann_join_tables --input {params.tmpdir} --output {output.abund} &>> {log}
    rm {params.tmpdir}/*
    
    for file in {input.cov}; do ln -s $(pwd)/$file $(pwd)/{params.tmpdir} &>> {log}; done
    humann_join_tables --input {params.tmpdir} --output {output.cov} &>> {log}
    rm {params.tmpdir}/*
    """

rule humann_map:
  name: "prok-annotate.smk HUMAnN3 map results"
  input:
    tsv=relpath("annotate/prok/output/genefamilies_merged.tsv"),
    ecdb=os.path.join(config['humann-db'], "utility_mapping/map_level4ec_uniref90.txt.gz"),
    eggnogdb=os.path.join(config['humann-db'], "utility_mapping/map_eggnog_uniref90.txt.gz"),
    godb=os.path.join(config['humann-db'], "utility_mapping/map_go_uniref90.txt.gz"),
    kodb=os.path.join(config['humann-db'], "utility_mapping/map_ko_uniref90.txt.gz"),
    pfamdb=os.path.join(config['humann-db'], "utility_mapping/map_pfam_uniref90.txt.gz"),
    unirefnamedb=os.path.join(config['humann-db'], "utility_mapping/map_uniref90_name.txt.bz2"),
    ecnamedb=os.path.join(config['humann-db'], "utility_mapping/map_ec_name.txt.gz"),
    eggnognamedb=os.path.join(config['humann-db'], "utility_mapping/map_eggnog_name.txt.gz"),
    gonamedb=os.path.join(config['humann-db'], "utility_mapping/map_go_name.txt.gz"),
    konamedb=os.path.join(config['humann-db'], "utility_mapping/map_ko_name.txt.gz"),
    pfamnamedb=os.path.join(config['humann-db'], "utility_mapping/map_pfam_name.txt.gz"),
  output:
    cpm=relpath("annotate/prok/output/genefamilies_merged-cpm.tsv"),
    rxn=relpath("annotate/prok/output/genefamilies_merged-cpm-metacycrxn.tsv"),
    ec=relpath("annotate/prok/output/genefamilies_merged-cpm-ec.tsv"),
    eggnog=relpath("annotate/prok/output/genefamilies_merged-cpm-eggnog.tsv"),
    go=relpath("annotate/prok/output/genefamilies_merged-cpm-go.tsv"),
    ko=relpath("annotate/prok/output/genefamilies_merged-cpm-ko.tsv"),
    pfam=relpath("annotate/prok/output/genefamilies_merged-cpm-pfam.tsv"),
    #tsvname=relpath("annotate/prok/output/genefamilies_merged-named.tsv"),
    rxnname=relpath("annotate/prok/output/genefamilies_merged-cpm-metacycrxn-named.tsv"),
    ecname=relpath("annotate/prok/output/genefamilies_merged-cpm-ec-named.tsv"),
    eggnogname=relpath("annotate/prok/output/genefamilies_merged-cpm-eggnog-named.tsv"),
    goname=relpath("annotate/prok/output/genefamilies_merged-cpm-go-named.tsv"),
    koname=relpath("annotate/prok/output/genefamilies_merged-cpm-ko-named.tsv"),
    pfamname=relpath("annotate/prok/output/genefamilies_merged-cpm-pfam-named.tsv"),
  params:
    outdir=relpath("annotate/prok/output"),
    basename=relpath("annotate/prok/output/genefamilies_merged-cpm"), 
    tmpdir=os.path.join(tmpd, "humann_map")
  conda: "../envs/biobakery3.yml"
  log: os.path.join(logdir, "humann_map.log")
  benchmark: os.path.join(benchmarks, "humann_map.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, threads, input: input.size_mb + 100000
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.outdir} {params.tmpdir}

    humann_renorm_table --input {input.tsv} --output {output.cpm} --units cpm --mode community --special y --update-snames &> {log}

    # humann_regroup_table --input {output.cpm} --output {output.rxn} --groups uniref90_rxn --function sum --ungrouped Y --protected Y &>> {log}
    humann_regroup_table --input {output.cpm} --output {output.ec} --custom {input.ecdb} --function sum --ungrouped Y --protected Y &>> {log}
    humann_regroup_table --input {output.cpm} --output {output.eggnog} --custom {input.eggnogdb} --function sum --ungrouped Y --protected Y &>> {log}
    humann_regroup_table --input {output.cpm} --output {output.go} --custom {input.godb} --function sum --ungrouped Y --protected Y &>> {log}
    humann_regroup_table --input {output.cpm} --output {output.ko} --custom {input.kodb} --function sum --ungrouped Y --protected Y &>> {log}
    humann_regroup_table --input {output.cpm} --output {output.pfam} --custom {input.pfamdb} --function sum --ungrouped Y --protected Y &>> {log}

    humann_rename_table --input {output.rxn} --output {output.rxnname} --names metacyc-rxn &>> {log}
    humann_rename_table --input {output.ec} --output {output.ecname} --custom {input.ecnamedb} &>> {log}
    humann_rename_table --input {output.eggnog} --output {output.eggnogname} --custom {input.eggnognamedb} &>> {log}
    humann_rename_table --input {output.go} --output {output.goname} --custom {input.gonamedb} &>> {log}
    humann_rename_table --input {output.ko} --output {output.koname} --custom {input.konamedb} &>> {log}
    humann_rename_table --input {output.pfam} --output {output.pfamname} --custom {input.pfamnamedb} &>> {log}
  
    rm -f {params.tmpdir}/*
    """


rule benchmark_summary:
  name: "prok-annotate.smk summaries performance benchmarks"
  input:
    os.path.join(benchmarks, "humann_merge.log")
  output:
    os.path.join(benchmarks, "summary.tsv")
  params:
    indir=benchmarks
  threads: 1
  shell:
    """
    cat {params.indir}/* | head -n 1 > {output}
    for file in $(find {params.indir}/*.log -type f); do echo -e "$(tail -n +2 $file)\t  $(basename $file)" >> {output}; done
    """
