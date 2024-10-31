import os

configdict = config['merge']


rule symlink_fastq:
  name: "merge.py symbolic link fastq"
