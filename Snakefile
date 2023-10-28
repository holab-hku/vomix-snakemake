import os
import glob
import sys
import platform
import pandas as pd

from source.utils.parse import parse_sample_list
from snakemake.utils import min_version, validate
from snakemake.exceptions import WorkflowError

workdir: config["workdir"]
configfile: "config.yml"

# Set temporary dir
if not os.getenv("TMPDIR"):
  os.environ["TMPDIR"] = "tmp"
  os.makedirs(os.environ["TMPDIR"],exist_ok=True)

# Set wildcard constraints
wildcard_constraints:
  sample_id = "[A-Za-z0-9_\-\.]+"

# Include rules
include: "pipeline/rules/viralcontigident.smk"

# Set output targets based on sample names

# 3) Viral Contig Identification
viralcontigident = expand("output/3__viralcontigident/{sample_id}/merged_scores.csv", sample_id = samples.keys())
