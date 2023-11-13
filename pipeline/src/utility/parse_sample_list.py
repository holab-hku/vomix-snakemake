#!/usr/bin/env python

import pandas as pd
import subprocess
import sys
import os

def validate_samples(samples):
    """
    Performs checks on the samples provided in the sample_list.tsv
    to either verify their existence locally or on the SRA.

    If you are using local files, you can either:
    1) Provide full file paths in the R1 and R2 columns of sample_list.tsv [only R1 if single-end]
    2) Place the fastq files in the config['datadir'] path with <sample>_{1,2}.fastq.gz naming format.
       If a file has no paired {2} label, it is assumed to be single-end.

    :param samples: samples dictionary
    :return:
    """
    for sample, items in samples.items():
        
        # Check if files exist locally
        if os.path.exists(samples[sample]["R1"]) and os.path.exists(samples[sample]["R2"]):
            continue

        # Check if it exists in SRA if not present locally
        acc = items['accession']
        cmd = ['sratools', 'info', acc]
        p = subprocess.run(cmd, stdout=subprocess.PIPE).stdout
        found = p.decode().split("\n")[0].split(":")[-1].lstrip()
        print(p)
        if int(found) == 0:
            sys.exit("""
########################## WARNING ###################################
# Accession: {} could not be found. Is it a valid SRA accession?     
#                                                                    #
# If you intend to use locally stored fastq files, make sure your    #
# sample file list contains columns named 'Read_file' and            #
# 'Paired_file' with paths pointing to the corresponding fastq files #
# for each sample.                                                   #
########################## WARNING ###################################
            """.format(acc))

def parse_sample_list(f, datadir):
    """
    Parse the sample list. Each sample is stored as a dictionary in the samples{} dictionary.
    samples{sample_name} will have the following information:

    samples[sample_name] = {'R1': 'path to R1',
                            'R2': 'path to R2',
                            'accession': 'accession id'}

    """

    samples = {}

    df = pd.read_csv(f, comment='#', header=0, sep='\t', index_col=0, dtype=str)
    df.fillna('', inplace=True)

    # If it's just a one-column file, expand it by assuming that the
    # SRA accessions are in the first column
    if df.shape[1] == 0 or 'accession' not in df.columns:
        df = df.assign(accession=pd.Series(df.index, index=df.index))
        df.index.name = "Sample"

    for sample in df.index:
        try:
            R1 = df.loc[sample, 'R1']
            R2 = df.loc[sample, 'R2']
        except KeyError:
            R1 = '{}_1.fastq.gz'.format(sample)
            R2 = '{}_2.fastq.gz'.format(sample)


        if 'accession' in df.columns:
            accession = df.loc[sample, 'accession']
        else:
            accession = ''
        
        samples[sample] = {'R1': '{dir}/{f}'.format(f=R1, dir=datadir),
                           'R2': '{dir}/{f}'.format(f=R2, dir=datadir),
                           'accession': accession}
    validate_samples(samples)
    
    return samples


# sample_list_file = './sample_list.txt'
# data_dir = '.'

# Call the parse_sample_list function to parse the sample list
#sample_dict = parse_sample_list(sample_list_file, data_dir)

