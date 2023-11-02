#!/bin/bash
#PBS -N TEST
#PBS -l nodes=1:ppn=12
#PBS -l mem=50g
#PBS -l walltime=03:00:00
#PBS -M eshekar@connect.hku.hk
# #PBS -m ae
#PBS -q cgsd
#PBS -o /home/eshekar/snakemake_pipelines/metaviromewrapper/out.qsub	
#PBS -e /home/eshekar/snakemake_pipelines/metaviromewrapper/err.qsub	



# RUN COMMANDS
cd /home/eshekar/snakemake_pipelines/metaviromewrapper 
source ~/.bashrc
conda activate metaviromewrapper
snakemake --use-conda --snakefile viralcontigident.smk output/3__viralcontigident/CRC_meta/merged_scores_filtered.csv -j 12 -F
