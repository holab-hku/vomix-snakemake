configfile: "config/preprocessing.yml"
logdir = relpath("preprocess/logs")
datadir = config["datadir"]
os.makedirs(logdir, exist_ok=True)


def retrieve_accessions(wildcards):
  try:
    acc=samples[wildcards.sample_id]['accession']
  except KeyError:
    acc=wildcards.sample_id
  return acc


rule download_fastq:
  name : "preprocessing.py download fastq from SRA"
  output:
    R1=os.path.join(datadir, "{sample_id}_1.fastq.gz"),
    R2=os.path.join(datadir, "{sample_id}_2.fastq.gz")
  params:
    download=config['dwnldparams'],
    pigz=config['pigzparams'],
    logdir=os.path.join(datadir, "/log"), 
    accessions= lambda wildcards: retrieve_accessions(wildcards),
    tmpdir="$TMPDIR/{sample_id}"
  log:
    os.path.join(datadir, "/log/{sample_id}.log")
  threads: 8
  shell:
    """
    mkdir -p {params.tmpdir}
    mkdir -p {params.logdir}

    fasterq-dump {params.accessions} \
        {params.download} \
        --split-3 \
        --skip-technical \
        --outdir {params.tmpdir} \
        --threads {threads} &> {log}

    pigz -p {threads} -c {params.pigz} {params.tmpdir}/*_1.fastq > {output.R1}
    pigz -p {threads} -c {params.pigz} {params.tmpdir}/*_2.fastq > {output.R2}

    rm -rf {params.tmpdir}

    """



rule fastp:
  name : "preprocessing.py fastp preprocess"
  input: 
    R1=os.path.join(datadir, '{sample_id}_1.fastq.gz'),
    R2=os.path.join(datadir, '{sample_id}_2.fastq.gz')
  output:
    R1=relpath("preprocess/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
    R2=relpath("preprocess/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz"), 
    html=relpath("preprocess/{sample_id}/report.fastp.html"),
    json=relpath("preprocess/{sample_id}/report.fastp.json")
  params:
    fastp=config['fastpparams'],
    tmpdir="$TMPDIR/{sample_id}"
  log: os.path.join(logdir, "fastp_{sample_id}.log")
  threads: 12
  conda: "../envs/fastp.yml"
  shell:
    """
    mkdir -p {params.tmpdir}

    fastp -i {input.R1} -I {input.R2} \
        -o {params.tmpdir}/R1.fastq.gz \
        -O {params.tmpdir}/R2.fastq.gz \
        --thread {threads} \
        --html {params.tmpdir}/tmp.html \
        --json {params.tmpdir}/tmp.json \
        {params.fastp} &> {log}

    mv {params.tmpdir}/R1.fastq.gz {output.R1}
    mv {params.tmpdir}/R2.fastq.gz {output.R2}
    mv {params.tmpdir}/tmp.html {output.html}
    mv {params.tmpdir}/tmp.json {output.json}

    rm -r {params.tmpdir}

    """


rule multiqc:
  name: "preprocessing.py preprocess report"
  input:
    logs=expand(relpath("preprocess/{sample_id}/report.fastp.json"), sample_id = samples.keys())
  output:
    relpath("reports/preprocess/preprocess_report.html"),
    relpath("reports/preprocess/preprocess_report_data/multiqc.log")
  params:
    searchdir=relpath("preprocess/"),
    outdir=relpath("reports/preprocess/"),
    tmpdir="$TMPDIR/multiqc"
  log: os.path.join(logdir, "multiqc.log")
  threads: 1
  conda: "../envs/multiqc.yml"
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    multiqc {params.searchdir} -f -o {params.tmpdir} -n preprocess_report.html 2> {log}
    mv {params.tmpdir}/preprocess_report* {params.outdir}

    rm -rf {params.tmpdir}
    """

