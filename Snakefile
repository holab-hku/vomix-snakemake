import os
import glob
import sys
import platform
import pandas as pd
import snakemake

#from snakemake import glob_wildcards
#from source.utils.parse import parse_sample_list
#from snakemake.utils import min_version, validate
#from snakemake.exceptions import WorkflowError

configfile: "config.yml"
workdir: config["workdir"]

### Set temporary dir
if not os.getenv("TMPDIR"):
  os.environ["TMPDIR"] = "tmp"
  os.makedirs(os.environ["TMPDIR"],exist_ok=True)

### Set wildcard constraints
wildcard_constraints:
  sample_id = "[A-Za-z0-9_\-\.]+"


### Parse sample list
df = pd.read_csv("sample_list.tsv", sep="\t", header=0)
sample_list = df.iloc[:, 0].tolist()


### Include rules
include: "viralcontigident.smk"


#### Set output targets based on sample names
viralcontigident = "output/3__viralcontigident/checkv/viruses.fna"



### Set rule all outputs

rule viralcontigident:
  input: viralcontigident



