import os 
import argparse
import pandas as pd




def extract_taxonomy(input_f, output_f):

	df = pd.read_csv(input_f, index_col=0)
	taxdf = df.loc[:, ['genomad_taxonomy']]
	taxdf['genomad_taxonomy'] = taxdf['genomad_taxonomy'].replace('Unclassified', 'Unassigned')

	# Assuming your DataFrame is called df
	taxdf[['virus', 'realm', 'kingdom', 'phylum', 'class', 'order', 'family']] = taxdf['genomad_taxonomy'].str.split(';', expand=True)
	taxdf = taxdf.drop('genomad_taxonomy', axis=1)

	taxdf.to_csv(output_f, index=True)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.description = "Takes the classification summary and the genomad taxonomy column and creates a taxonomy table"

	parser.add_argument("--input", dest="infile", help="input classification summary file", required=True)
	parser.add_argument("--output", dest="outfile", help="ouput csv file path", required=True)
	
	args = parser.parse_args()
	extract_taxonomy(args.infile, args.outfile)



