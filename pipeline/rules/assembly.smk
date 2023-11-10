#localrules:
#all

###########
# MEGAHIT #
###########

rule megahit:
  input:
    R1 = "output/1__preprocessing/{sample_id}_R1.cut.trim.fastq.gz",
    R2 = "output/1__preprocessing/{sample_id}_R2.cut.trim.fastq.gz"
  output:
    fasta = "output/2__assembly/{sample_id}/final.contigs.fa",
    R1 =  "output/1__preprocessing/{sample_id.fastq.gz", 
    R2 = "output/1__preprocessing/{sample_id}.fastq.gz"
  params:
    megahitparams = config['megahitparams'],
    min_contig_len = config["megahit_min_contig_len"],
    out_dir = "output/2__assembly/{sample_id}",
    tmp_dir = "$TMPDIR/{sample_id}",
  log: "logs/2__assembly_cdhitderep.log"
  threads: 10
  shell:
    """
    wd=$(pwd)
    ln -s $wd/{input.R1} $wd/{output.R1}
    ln -s $wd/{input.R2} $wd/{output.R2}

    mkdir -p {params.tmp_dir}
    megahit -1 {input.R1} -2 {input.R2} \
        --min-contig-len {params.min_contig_len} -o {params.tmp_dir} \
        -t {threads} {params.megahitparam} > {log} 2>&1
    mv {params.tmp_dir}/final.contigs.fa {output.fa}
    mv {params.tmp_dir}/opts.txt {params.out_dir}
    
    """


###################
# MetaviralSPAdes #
###################


