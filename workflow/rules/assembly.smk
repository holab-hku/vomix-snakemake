localrules:
  symlink_assembly

###########
# MEGAHIT #
###########
rule symlink_assembly:
  input:
    R1 = "results/preprocess/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz",
    R2 = "results/preprocess/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz"
  output:
    R1 =  "results/assembly/{sample_id}/output/{sample_id}_R1.fastq.gz",
    R2 = "results/assembly/{sample_id}/output/{sample_id}_R2.fastq.gz"
  params:
    output_dir = "results/assembly/{sample_id}/output"
  shell:
    """
    wd=$(pwd)
    ln -s $wd/{input.R1} $wd/{output.R1}
    ln -s $wd/{input.R2} $wd/{output.R2}
    """


rule megahit:
  input:
    R1 = "results/assembly/{sample_id}/output/{sample_id}_R1.fastq.gz",
    R2 = "results/assembly/{sample_id}/output/{sample_id}_R2.fastq.gz"
  output:
    fasta = "results/assembly/{sample_id}/output/final.contigs.fa"
  params:
    megahitparams = config['megahitparams'],
    min_contig_len = config["megahit_min_contig_len"],
    inter_dir = "results/assembly/{sample_id}/intermediate/megahit",
    tmp_dir = "$TMPDIR/{sample_id}"
  log: "logs/assembly_{sample_id}_megahit.log"
  conda: "../envs/megahit.yml"
  threads: 48
  shell:
    """
    rm -rf {params.tmp_dir} {params.inter_dir} {output}
    mkdir -p {params.inter_dir}

    megahit -1 {input.R1} -2 {input.R2} \
        --min-contig-len {params.min_contig_len} \
        -o {params.tmp_dir} \
        -t {threads} \
        {params.megahitparams} &> {log} 

    mv {params.tmp_dir}/final.contigs.fa {output.fasta}
    mv {params.tmp_dir}/* {params.inter_dir}
    
    """


rule assembly_stats:
  input:
    expand("results/assembly/{sample_id}/output/final.contigs.fa", sample_id = samples.keys())
  output:
    stats = "results/report/assembly/assemblystats.tsv",
    sizedist = "results/report/assembly/assembly_size_dist.tsv"
  params:
    output_dir = "results/report/assembly",
    names = "test"
  log: "logs/assembly_assemblystats.log"
  conda: "../envs/utility.yml"
  shell:
    """
    mkdir -p {params.output_dir}
    python pipeline/src/assembly/assemblystats.py -i {input} -n {params.names} --size-dist-file {output.sizedist} > {output.stats}
     """


###################
# MetaviralSPAdes #
###################

 





