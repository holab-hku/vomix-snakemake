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

samples, assemblies = parse_sample_list(config["samplelist"], config['datadir'])
accession_ids = [sample['accession'] for sample in samples.values()]


### Include rules
include: "preprocessing.smk"
include: "viralcontigident.smk"


#### Set output targets based on sample names
preprocess = expand("output/preprocess/{sample_id}/output/{sample_id}_R{i}_cut.trim.filt.fastq.gz", sample_id = samples.keys(), i = [1,2])
preprocess += ["output/report/preprocess/preprocess_report.html"]

viralcontigident = expand("output/viralcontigident/{sample_id}/output/viral.contigs.fa", sample_id = samples.keys())
viralcontigident += "output/viralcontigident/output/checkv/viruses.fna"




### Set rule all outputs
rule preprocess:
  input: preprocess
rule viralcontigident:
  input: viralcontigident



