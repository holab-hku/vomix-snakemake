rule blastp_taxonomy:
  input: "
  output:
  params:
  conda: "../envs/blast.yml"
  log:
  threads:
  shell:
    """
    """


rule vcontact2_taxonomy:
  input:
  output:
  params:
  conda: "../envs/vcontact2.yml"
  log:
  threads:
  shell:
    """
    """

rule viphogs_taxonomy:
  input:
  output:
  params:
  conda: "../envs/viphogs.yml"
  log:
  threads:
  shell:
    """
    """

rule phagcn_taxonomy:
  input:
  output:
  params:
  conda: "../envs/phabox.yml"
  log:
  threads:
  shell:
    """
    """

rule merge_taxonomy:
  input:
  output:
  params:
  conda:
  log:
  threads:
  shell:
    """
    """

rule consensus_taxonomy:
  input:
  output:
  params:
  conda:
  log:
  threads:
  shell:
    """
    """
