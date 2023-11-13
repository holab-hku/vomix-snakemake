import os
import glob
import sys
import platform
import pandas as pd
import snakemake

#from snakemake import glob_wildcards
from pipeline.src.utility.parse_sample_list import parse_sample_list
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


### Parse sample python dictionary

### It has the format samples[sample_name]  =  {'R1': 'path to R1',
#                                               'R2': 'path to R2',
#                                               'accession': 'accession id'}

samples = parse_sample_list(config["samplelist"], config['datadir'])
accession_ids = [sample['accession'] for sample in samples.values()]


### Include rules
include: "preprocessing.smk"
include: "viralcontigident.smk"


#### Set output targets based on sample names
# preprocess = expand(os.path.join(config['datadir'], "{sample_id}_{i}.fastq.gz"), sample_id = accession_ids, i = [1,2])
preprocess = expand("output/preprocess/{sample_id}_R{i}.cut.trim.fastq.gz", sample_id = samples.keys(), i = [1,2])
viralcontigident = "output/3__viralcontigident/checkv/viruses.fna"



### Set rule all outputs
#rule preprocess:
#  input: preprocess
rule viralcontigident:
  input: viralcontigident



