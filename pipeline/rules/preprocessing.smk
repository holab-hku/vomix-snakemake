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

