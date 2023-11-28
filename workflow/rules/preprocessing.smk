

# a lambda function to get sample accesions 
def retrieve_accessions(wildcards):
  try:
    acc = samples[wildcards.sample_id]['accession']
  except KeyError:
    acc = wildcards.sample_id
  print(acc)
  return acc


rule download_fastq:
  output:
    R1 = os.path.join(config["datadir"], "{sample_id}_1.fastq.gz"),
    R2 = os.path.join(config["datadir"], "{sample_id}_2.fastq.gz")
  params:
    dwnldparams = config['dwnldparams'],
    pigzparams = config['pigzparams'],
    log_dir = os.path.join(config["datadir"], "/log"), 
    accessions= lambda wildcards: retrieve_accessions(wildcards),
    tmp_dir = "$TMPDIR/{sample_id}"
  log:
    os.path.join(config["datadir"], "/log/{sample_id}.log")
  threads: 8
  shell:
    """
    mkdir -p {params.tmp_dir}
    mkdir -p {params.log_dir}

    fasterq-dump {params.accessions} \
        {params.dwnldparams} \
        --split-3 \
        --skip-technical \
        --outdir {params.tmp_dir} \
        --threads {threads} &> {log}

    pigz -p {threads} -c {params.pigzparams} {params.tmp_dir}/*_1.fastq > {output.R1}
    pigz -p {threads} -c {params.pigzparams} {params.tmp_dir}/*_2.fastq > {output.R2}

    rm -rf {params.tmp_dir}

    """



rule fastp:
  input: 
    R1 = os.path.join(config['datadir'], '{sample_id}_1.fastq.gz'),
    R2 = os.path.join(config['datadir'], '{sample_id}_2.fastq.gz')
  output:
    R1 = "results/preprocess/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz",
    R2 = "results/preprocess/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz", 
    html = "results/preprocess/{sample_id}/report.fastp.html", 
    json = "results/preprocess/{sample_id}/report.fastp.json"
  params:
    fastpparams = config['fastpparams'],
    R1_tmp = "$TMPDIR/{sample_id}/R1.fastq.gz",
    R2_tmp = "$TMPDIR/{sample_id}/R2.fastq.gz",
    html_tmp = "$TMPDIR/{sample_id}/fastp.html",
    json_tmp = "$TMPDIR/{sample_id}/fastp.json",
    tmp_dir = "$TMPDIR/{sample_id}/"
  log: "logs/preprocess_{sample_id}_fastp.log"
  threads: 12
  conda: "../envs/fastp.yml"
  shell:
    """
    mkdir -p {params.tmp_dir}

    fastp -i {input.R1} -I {input.R2} \
        -o {params.R1_tmp} -O {params.R2_tmp} \
        --thread {threads} \
        --html {params.html_tmp} \
        --json {params.json_tmp} \
        {params.fastpparams} &> {log}

    mv {params.R1_tmp} {output.R1}
    mv {params.R2_tmp} {output.R2}
    mv {params.html_tmp} {output.html}
    mv {params.json_tmp} {output.json}

    rm -r {params.tmp_dir}

    """


rule multiqc:
  input:
    fastplogs = expand("results/preprocess/{sample_id}/report.fastp.json", sample_id = samples.keys())
  output:
    "workflow/report/preprocess/preprocess_report.html",
    "workflow/report/preprocess/preprocess_report_data/multiqc.log"
  params:
    search_dir = "results/preprocess/",
    output_dir = "workflow/report/preprocess/",
    tmp_dir = "$TMPDIR/multiqc/"
  log: "logs/preprocess_multiqc.log"
  threads: 1
  conda: "../envs/multiqc.yml"
  shell:
    """
    rm -rf {params.tmp_dir}
    rm -rf {params.output_dir}
    mkdir -p {params.tmp_dir}
    mkdir -p {params.output_dir}

    multiqc {params.search_dir} -f -o {params.tmp_dir} -n preprocess_report.html 

    mv {params.tmp_dir}preprocess_report* {params.output_dir}

    """

