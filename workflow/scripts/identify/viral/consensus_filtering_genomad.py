import os 
import argparse
import pandas as pd 


def consensus_filtering(classification_summary, genomad_min, summary_out, provirus_out, virus_out):

	df = pd.read_csv(classification_summary, index_col=0)

	# add a new column called hit-score which shows how many softwares it has a hit from 
	# according to the config file thresholds set

	df['hit-score'] = 0
	df.loc[(df['genomad_virus_score'] > genomad_min), 'hit-score'] += 1
	
	df['provirus-hit-score'] = 0
	df.loc[(df['checkv_provirus'] == 'Yes'), 'provirus-hit-score'] += 1
	df.loc[(df['genomad_topology'] == 'Provirus'), 'provirus-hit-score'] += 1
	
	# set sequences identified as proviruses by both checkV and genoad as true pro-viruses
	# if only one tool says it's a provirus, it will label as lower confidence
	df['type'] = 'Virus'
	df.loc[(df['provirus-hit-score'] == 2), 'type'] = 'Provirus'
	df.loc[(df['provirus-hit-score'] == 1), 'type'] = 'Provirus (Low-Confidence)'

	# filter the dataframe if not 'type' is 'Virus' but contamination is higher than 10%
	df = df[(df['type'] != 'Virus') | (df['checkv_contamination'] <= 10)]
	
	# do rigorous filtering based on whether checkv completeness  assessment 
	# The filtering is done through a loop for code readability reasons
	filtered_rows = []
	
	for index, row in df.iterrows():
		# 1) Keep every sequence where the  CheckV genome is set as Complete or High-quality
		if row['checkv_quality'] in ('High-quality', 'Complete'):
			filtered_rows.append(row)

		# 4) If genome is medium or low-quality or not-determined, keep if genomad score is very high
		# (according to secondary geNomad threshold) 
		else: 
			if row['hit-score'] == 1:
				filtered_rows.append(row)
			else:
				continue
			
	filtered_df = pd.concat(filtered_rows, axis=1).transpose()

	# drop any sequence where more than one proviruse region has been identified by geNomad
	# but checkV does not see a provirus
	#filtered_df = filtered_df[~(filtered_df.index.duplicated()) & (filtered_df['provirus-hit-score'] != 2)]
	
	# drop duplicate indices 
	# these regions represent the ones where geNomad found two proviruses 
	# we can find a cleaner way to handle this later
	filtered_df = filtered_df[~filtered_df.index.duplicated(keep='first')]

	# make all indecies unique 
	provirusdf = filtered_df[(filtered_df['type'] == 'Provirus') | (filtered_df['type'] == 'Provirus (Low-Confidence)')]
	virusdf = filtered_df[filtered_df['type'] == 'Virus']


	# save outputs
	filtered_df.to_csv(summary_out, index=True)
	provirusdf.index.to_frame().to_csv(provirus_out, header=False, index=False, sep='\t')
	virusdf.index.to_frame().to_csv(virus_out, header=False, index=False, sep='\t')


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Perform consensus viral classification on the outputs to find vOTUs')

	parser.add_argument('--classification_results', type=str, help='Path to final combined classification summary')
	parser.add_argument('--genomad_min_score', type=float, help='Minimum threshold for geNomad to be considered as a hit')
	parser.add_argument('--summary_out', type=str, help='Output file path to final filtered summary of vOTUs')
	parser.add_argument('--provirus_list', type=str, help='Output file path to be list of proviral sequences')
	parser.add_argument('--virus_list', type=str, help='Output file path to be list of viral (non-proviral) sequences')

	args = parser.parse_args()

	consensus_filtering(args.classification_results, args.genomad_min_score, 
			args.summary_out, args.provirus_list, args.virus_list)

