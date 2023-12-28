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

		# check if files exist locally
		if os.path.exists(samples[sample]["R1"]) and os.path.exists(samples[sample]["R2"]):
			continue

		# check if it exists in SRA if not present locally
		acc = items['accession']
		cmd = ['efetch', '-db', 'sra', '-id', acc, '-format', 'runinfo']
		#cmd = ['sratools', 'info', acc]
		#p = subprocess.run(cmd, stdout=subprocess.PIPE).stdout
		#found = p.decode().split("\n")[0].split(":")[-1].lstrip()
		found = 1
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
	
	df = pd.read_csv(f, comment='#', header=0, sep='\t', index_col=None, dtype=str)
	df = df.replace(r'^\s*$', float('nan'), regex=True)
	if 'assembly' not in df.columns:
		    df['assembly'] = float('nan')

	# iterate through df and if sample_id is missing, replace is with SRA accession
	# also add assembly as sample_id if no assembly is given (single-sample assembly)
	for i, row in df.iterrows():
		if pd.isna(row['sample_id']) and not pd.isna(row['accession']):
			df.loc[i, 'sample_id'] = row['accession']

		elif pd.isna(row['sample_id']) and pd.isna(row['accession']):
			raise ValueError("Column with both empty sample_id and SRA accession, please provide at least one")
	
	# if there are no assemblies set, make it the same as sample_id (single sample assembly)
	df['assembly'] = [row['assembly'] if not pd.isna(row['assembly']) else row['sample_id'] for _, row in df.iterrows()]
	
	# if there are no R1 or R2s, generate them
	df['R1'] = [row['R1'] if not pd.isna(row['R1']) else '{dir}{s}_1.fastq.gz'.format(dir= datadir, s=row['sample_id']) for _, row in df.iterrows()]
	df['R2'] = [row['R2'] if not pd.isna(row['R2']) else '{dir}{s}_2.fastq.gz'.format(dir = datadir, s=row['sample_id']) for _, row in df.iterrows()]
		
	df.fillna('', inplace=True)
	
	# set unique names for the file index
	try:
		df.set_index('sample_id', inplace=True)
	except ValueError:
		print("sample_id error. Values in the sample_id column may be non-unique")


	

	# If it's just a one-column file, expand it by assuming that the
	# SRA accessions are in the first column
	if (df.shape[1] == 0):
		if ('accession' not in df.columns):
			df = df.assign(accession = pd.Series(df.index, index = df.index))
			df.index.name = "sample_id"
	
	# Remove duplicates rows and throw a warning
	duplicates = df[df.duplicated()]
	if not duplicates.empty:
		print("Duplicate samples found in sample_list.tsv. Automatically deleting them. Please beware of duplicate samples and typos.")
		df.drop_duplicates(inplace=True)


	# setup co-assemblies if needed
	assemblies = {}
	for assembly in df.assembly.unique():
		df_filt = df[df['assembly'] == assembly]
		R1s = df_filt['R1'].tolist()
		R2s = df_filt['R2'].tolist()

		try:
			accessions = df_filt['accession'].tolist()
		except (NameError, AttributeError, KeyError):
			accessions = []

		sample_ids = df_filt.index.tolist()

		assemblies[assembly] = {'R1': R1s, 'R2': R2s, 
				'sample_id' : sample_ids, 
				'accession' : accessions}

	for sample_id in df.index:
		try:
			R1 = df.loc[sample_id, 'R1']
			R2 = df.loc[sample_id, 'R2']
		except KeyError:
			R1 = '{}_1.fastq.gz'.format(sample_id)
			R1 = '{dir}/{f}'.format(f=R1, dir=datadir)
			R2 = '{}_2.fastq.gz'.format(sample_id)
			R2 = '{dir}/{f}'.format(f=R2, dir=datadir)
 		
		if 'accession' in df.columns:
			accession = df.loc[sample_id, 'accession']
		else:
			accession = ''

		samples[sample_id] = {'R1': R1,
				'R2': R2,
				'accession': accession}

	
	validate_samples(samples)
	return samples, assemblies
