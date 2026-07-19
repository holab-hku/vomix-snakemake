import math

logdir = relpath("identify/viral/logs")
tmpd = relpath("identify/viral/tmp")
benchmarks=relpath("identify/viral/benchmarks")

os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)

clustering_iter = config["cluster-iter"] # nLayers
n_chunks_layer_1 =  2 ** (clustering_iter - 1) # nCluster Chunks for Layer 1

#### Helper function to choose inputs for recursive clustering
def get_cdhit_inputs(wildcards):

    """Calculates which files to merge based on the current layer and chunk."""
    layer = int(wildcards.layer)
    chunk = int(wildcards.chunk)
    
    if layer == 1:
        # Layer 1 just grabs the raw split pieces
        return f"{tmpd}/splits/chunk_{chunk}.fa"
    else:
        # Layer > 1 grabs the two specific clustered chunks from the previous layer
        prev_layer = layer - 1
        child1 = chunk * 2
        child2 = chunk * 2 + 1
        return [f"{tmpd}/layer_{prev_layer}/chunk_{child1}.fa", f"{tmpd}/layer_{prev_layer}/chunk_{child2}.fa"]
        
    




### Read single fasta file if input
if config['fasta'] != "" and config["module"] == "cluster-fast":
  fastap = readfasta(config["fasta"])
  sample_id = config["sample-name"]
  assembly_ids = [sample_id]
else:
  fastap = relpath("identify/viral/intermediate/scores/combined.viralcontigs.fa")
  sample_id = "combined.viralcontigs"
  assembly_ids = [sample_id]

### MASTER RULE
rule done_log:
  name: "clustering.smk Done. removing tmp files"
  localrule: True
  input:
    relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa.clstr"), 
  output:
    os.path.join(logdir, "clustering-done.log")
  shell: "touch {output}"


### RULES 
# 1) RAPID CLUSTERING (clustering-fast=True)
if config["clustering-fast"]:
  rule makeblastdb_derep:
    name: "clustering.smk make blast db [--clustering-fast]"
    input: 
      fastap
    output: 
      expand(relpath("identify/viral/intermediate/derep/db.{suffix}"), suffix=["ntf", "ndb"])
    params:
      outdir=relpath("identify/viral/intermediate/derep/"), 
      dbtype='nucl', 
      tmpdir=tmpd
    log: os.path.join(logdir, "clustering/makeblastdb.log")
    benchmark: os.path.join(benchmarks, "makeblastdb.log")
    conda: "../envs/checkv.yml"
    threads: 1
    resources:
      mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
    shell:
      """
      rm -rf {params.tmpdir}/* {params.outdir}
      mkdir -p {params.tmpdir} {params.outdir}

      makeblastdb -in {input} -dbtype {params.dbtype} -out {params.tmpdir}/db &> {log}

      mv {params.tmpdir}/* {params.outdir}
      rm -rf {params.tmpdir}/*

      """

  rule megablast_derep:
    name: "clustering.smk megablast [--clustering-fast]"
    input:
      fasta=fastap, 
      dbcheckpoints=expand(relpath("identify/viral/intermediate/derep/db.{suffix}"), suffix=["ntf", "ndb"])
    output:
      relpath("identify/viral/intermediate/derep/blast_out.csv")
    params:
      db=relpath("identify/viral/intermediate/derep/db"),
      outfmt="'6 std qlen slen'",
      maxtargetseqs=10000, 
      tmpdir=tmpd
    log: os.path.join(logdir, "megablastpairwise.log")
    benchmark: os.path.join(benchmarks, "megablastpairwise.log")
    conda: "../envs/checkv.yml"
    threads: 64
    resources:
      mem_mb=lambda wildcards, attempt: attempt * 72 * 10**3
    shell: 
      """
      rm -rf {params.tmpdir}/*
      mkdir -p {params.tmpdir}
      
      blastn -query {input.fasta} \
          -db {params.db} \
          -outfmt {params.outfmt} \
          -max_target_seqs {params.maxtargetseqs} \
          -out {params.tmpdir}/tmp.csv \
          -num_threads {threads} &> {log}
    
      mv {params.tmpdir}/tmp.csv {output}
      rm -rf {params.tmpdir}/*

      """
  
  rule anicalc_derep:
    name : "clustering.smk calculate ani [--clustering-fast]"
    input:
      relpath("identify/viral/intermediate/derep/blast_out.csv")
    output: 
      relpath("identify/viral/intermediate/derep/ani.tsv")
    params:
      script="workflow/scripts/clust_anicalc.py", 
      tmpdir=tmpd
    log: os.path.join(logdir, "anicalc.log")
    benchmark: os.path.join(benchmarks, "anicalc.log")
    conda: "../envs/checkv.yml"
    threads: 1
    resources:
      mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
    shell:
      """
      rm -rf {params.tmpdir}/* 
      mkdir -p {params.tmpdir}

      python {params.script} \
          -i {input} \
          -o {params.tmpdir}/tmp.tsv &> {log}

      mv {params.tmpdir}/tmp.tsv {output}
      rm -rf {params.tmpdir}/*
      """


  rule aniclust_derep:
    name : "clustering.smk cluster [--clustering-fast]"
    input:
      fa=fastap, 
      ani=relpath("identify/viral/intermediate/derep/ani.tsv")
    output:
      tsv= relpath("identify/viral/output/derep/clusters.tsv"),
      reps=relpath("identify/viral/output/derep/cluster_representatives.txt")
    params:
      script="workflow/scripts/clust_ani.py",
      minani=config["vOTU-ani"],
      targetcov=config["vOTU-targetcov"],
      querycov =config["vOTU-querycov"], 
      tmpdir=tmpd
    log: os.path.join(logdir, "aniclust.log")
    benchmark: os.path.join(benchmarks, "aniclust.log")
    conda: "../envs/checkv.yml"
    threads: 1
    resources:
      mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
    shell:
      """
      rm -rf {params.tmpdir}/*
      mkdir -p {params.tmpdir}

      python {params.script} \
          --fna {input.fa} \
          --ani {input.ani} \
          --out {params.tmpdir}/tmp.tsv \
          --min_ani {params.minani} \
          --min_tcov {params.targetcov} \
          --min_qcov {params.querycov} &> {log}
      mv {params.tmpdir}/tmp.tsv {output.tsv}
      cut -f1 {output.tsv} > {output.reps}

      rm -rf {params.tmpdir}/*
      """


  rule filtercontigs_derep:
    name: "clustering.smk filter dereplicated viral contigs"
    input: 
      fna=fastap, 
      reps=relpath("identify/viral/output/derep/cluster_representatives.txt")
    output:
      relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa")
    params:
      outdir=relpath("identify/viral/checkv/output"),
      tmpdir=tmpd
    log: os.path.join(logdir, "filterderep.log")
    conda: "../envs/seqkit-biopython.yml" 
    threads: 1
    resources:
      mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
    shell:
      """
      rm -rf {params.tmpdir}/*
      mkdir -p {params.tmpdir}

      seqkit grep {input.fna} -f {input.reps} > {params.tmpdir}/tmp.fa 2> {log}
      mv {params.tmpdir}/tmp.fa {output}

      rm -rf {params.tmpdir}/*
      """


# 2) CD-HIT CLUSTERING (clustering-fast=False)
else:
  rule cdhit_split_input:
    name: "clustering.smk Split input fasta"
    input:
        fastap
    output:
        expand(os.path.join(tmpd, "splits", f"chunk_{{chunk}}.fa"), chunk=range(n_chunks_layer_1))
    params:
        pieces = n_chunks_layer_1,
        outdir = os.path.join(tmpd, "splits")
    log: os.path.join(logdir, "clustering/split_input.log")
    benchmark: os.path.join(benchmarks, "split_input.log")
    conda: "../envs/seqkit-biopython.yml"
    threads: 1
    shell:
      """
        rm -rf {params.outdir}/*
        mkdir -p {params.outdir}
        
        seqkit split2 {input} -p {params.pieces} -O {params.outdir}/
        
        # Converts seqkit's default padded names (e.g., .part_001.fa) 
        # to our strict 0-indexed names (chunk_0.fa, chunk_1.fa) so math works.
        counter=0

        # Safely collect the files without triggering Snakemake or Bash syntax errors
        set -- {params.outdir}/*.fa {params.outdir}/*.fna {params.outdir}/*.fasta
        for file in "$@"; do
          ext="${{file##*.}}"
          mv "$file" "{params.outdir}/chunk_${{counter}}.${{ext}}"
          counter=$((counter+1))
        done        
        """
        
  rule cdhit_recursive_cluster:
    name: "clustering.smk CD-Hit recursive clustering [clustering-fast=False]"
    input:
        get_cdhit_inputs
    output:
        fa = os.path.join(tmpd, f"layer_{{layer}}", f"chunk_{{chunk}}.fa"),
        clstr = os.path.join(tmpd, f"layer_{{layer}}", f"chunk_{{chunk}}.fa.clstr")
    params:
        cdhitparams = config['cdhit-params']
    threads: 32
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 72 * 10**3
    conda: "../envs/cd-hit.yml"
    log: os.path.join(logdir, "clustering/cdhit_layer_{layer}_chunk_{chunk}.log")
    benchmark: os.path.join(benchmarks, "cdhit_layer_{layer}_chunk_{chunk}.log")
    shell:
        """
        mkdir -p $(dirname {output.fa})
        
        if [ "{wildcards.layer}" -eq "1" ]; then
            # LAYER 1: {input} contains exactly one file
            cd-hit -i {input} -o {output.fa} -T {threads} {params.cdhitparams} &> {log}
            
        else
            # LAYER > 1: {input} automatically expands to "fileA fileB" 
            # 'cat' will naturally pool them both into the temp file
            cat {input} > {output.fa}.tmp_input
            cd-hit -i {output.fa}.tmp_input -o {output.fa} -T {threads} {params.cdhitparams} &> {log}
            # Clean up the pooled temporary file to save disk space
            rm {output.fa}.tmp_input
        fi
        """

  rule cdhit_derep_finalize:
    name: "clustering.smk Create final dereplicated dataset"
    input:
        # We explicitly request Chunk 0 of the Final Layer (the bottom of the tree)
        fa = os.path.join(tmpd, f"layer_{clustering_iter}", "chunk_0.fa"),
        clstr = os.path.join(tmpd, f"layer_{clustering_iter}", "chunk_0.fa.clstr")
    output:
        fa = relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa"),
        clstr = relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa.clstr"), 
        done = os.path.join(logdir, "clustering-sensitive-done.log")
    log: os.path.join(logdir, "clustering/viral_cdhit.log")
    benchmark: 
        os.path.join(benchmarks, "identify/viral_cdhit.log")
    shell:
        """
        cp {input.fa} {output.fa}
        cp {input.clstr} {output.clstr}
        
        touch {output.done}
        """