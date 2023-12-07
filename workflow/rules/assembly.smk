configfile: "config/assembly.yml"

rule symlink_assembly:
  name : "assembly.py create symbolic links"
  localrule: True
  input:
    R1="results/preprocess/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz",
    R2="results/preprocess/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz"
  output:
    R1= "results/assembly/{sample_id}/output/{sample_id}_R1.fastq.gz",
    R2="results/assembly/{sample_id}/output/{sample_id}_R2.fastq.gz"
  params:
    outdir="results/assembly/{sample_id}/output"
  shell:
    """
    mkdir -p {params.outdir}
    wd=$(pwd)

    ln -s $wd/{input.R1} $wd/{output.R1}
    ln -s $wd/{input.R2} $wd/{output.R2}

    """

rule megahit:
  name : "assembly.py MEGAHIT assembly"
  input:
    R1="results/assembly/{sample_id}/output/{sample_id}_R1.fastq.gz",
    R2="results/assembly/{sample_id}/output/{sample_id}_R2.fastq.gz"
  output:
    fasta="results/assembly/{sample_id}/output/final.contigs.fa"
  params:
    parameters=config['megahitparams'],
    minlen=config["megahit_min_contig_len"],
    outdir="results/assembly/{sample_id}/output/",
    interdir="results/assembly/{sample_id}/intermediate/megahit",
    tmpdir="$TMPDIR/{sample_id}"
  log: "logs/assembly_{sample_id}_megahit.log"
  conda: "../envs/megahit.yml"
  threads: 48
  shell:
    """
    rm -rf {params.tmpdir} {params.interdir} {params.outdir}/*
    mkdir -p {params.interdir} {params.outdir}

    megahit -1 {input.R1} -2 {input.R2} \
        --min-contig-len {params.minlen} \
        -o {params.tmpdir} \
        -t {threads} \
        {params.parameters} &> {log} 

    mv {params.tmpdir}/final.contigs.fa {output.fasta}
    mv {params.tmpdir}/* {params.interdir}
    
    """


rule assembly_stats:
  name: "assembly.py aggregate assembly statistics"
  input:
    expand("results/assembly/{sample_id}/output/final.contigs.fa", sample_id=samples.keys())
  output:
    stats="workflow/report/assembly/assemblystats.tsv",
    sizedist="workflow/report/assembly/assembly_size_dist.tsv"
  params:
    script="workflow/scripts/assembly/assemblystats.py",
    outdir="workflow/report/assembly", 
    tmpdir="$TMPDIR/{sample_id}"
  log: "logs/assembly_assemblystats.log"
  conda: "../envs/utility.yml"
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}/* 
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} -i {input} --size-dist-file {params.tmpdir}/tmp1.tsv > {params.tmpdir}/tmp2.tsv 2> {log} 
    
    mv {params.tmpdir}/tmp1.tsv {output.sizedist}
    mv {params.tmpdir}/tmp2.tsv {output.stats}

     """

