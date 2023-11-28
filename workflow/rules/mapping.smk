rule bowtie_build:
  input:
    "output/assembly/{sample_id}/output/final.contigs.fa"
  output:
    expand("output/mapping/{sample_id}/output/final.contigs.fa.{index}.bt2l", index = range(1,5))
  params:
    prefix = "output/mapping/{sample_id}/output/final.contigs.fa",
    tmp_dir = "$TMPDIR/{sample_id}"
  log: = "logs/assembly_{sample_id}_bowtiebuild.log"
  conda: "pipeline/envs/bowtie2.yml"
  shell:
    """
    bowtie2-build --large-index --threads {threads} {input} {params.prefix}
    """


rule symlink_bowtie_index:
  input: 
    expand("output/mapping/{sample_id}/output/final.contigs.fa.{index}.bt2l", index = range(1,5))
  output:
    expand("output/assembly/{sample_id}/output/final.contigs.fa.{index}.bt2l", index = range(1,5))
  params:
    output_dir = "output/assembly/{sample_id}/output/"
  shell: "ln -s {input} {params.output_dir}"



rule bowtie2:
  input:
    R1= "output/assembly/{sample_id}/output/{sample_id}_R1.fastq.gz",
    R2 = "output/assembly/{sample_id}/output/{sample_id}_R2.fastq.gz",
    bowtie = expand("output/assembly/{sample_id}/output/final.contigs.fa.{index}.bt2l", index = range(1,5))
  output:
    "results/mapping/{sample_id}/output/{sample_id}.bam"
  params:
    bowtie2params = config["bowtie2params"],
    prefix = "results/mapping/{sample_id}/final.contigs.fa",
    tmp_dir = "$TMPDIR/{sample_id}",
    tmp_out = "$TMPDIR/{sample_id}/{sample_id}.bam"
  log: "logs/mapping_{sample_id}_bowtie2.log"
  conda: "pipeline/envs/bowtie2.yml"
  threads: 16
  shell:
    """
    rm -rf {params.tmp_dir}
    mkdir -p {params.tmp_dir}

    bowtie2 \
        -1 {input.R1} -2 {input.R2} \
        -p {threads} \
        -x {params.prefix} \
        {params.bowtie2params}  2>{log} | samtools view -h -b -o | samtools sort - -o {params.tmp_out} 2> {log}
    
    mv {params.tmp_out} {output}
    """


rule multiqc_bowtie2:
  input:
    bowtie2logs = expand("logs/mapping_{sample_id}_bowtie2.log", sample_id = samples.keys()),
    bam = expand("results/mapping/{sample_id}/output/{sample_id}.bam", sample_id = samples.keys())
  output:
    "output/report/mapping/preprocess_report.html",
    "output/report/mapping/preprocess_report_data/multiqc.log"
  params:
    search_dir = "output/logs/",
    output_dir = "output/report/mapping/",
    tmp_dir = "$TMPDIR/multiqc/"
  log: "logs/mapping_multiqc.log"
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
