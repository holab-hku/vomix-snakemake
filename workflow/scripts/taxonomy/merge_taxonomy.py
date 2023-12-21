import pandas as pd
from Bio import SeqIO
import numpy as np

def merge_taxonomy(diamond_out, viphogs_out, phagcn_out, genomad_out, contigs_fa, out_csv):
	""" Takes the output of different taxonomic annotation tools and returns a merged dataframe in csv format with column prefixes."""
	diamondout = diamond_out
	viphogsout = viphogs_out
	phagcnout = phagcn_out
	genomadout = genomad_out
	
	contigsfa = contigs_fa
	
	# Function
	headers = []
	with open(contigsfa, "r") as f:
		for record in SeqIO.parse(f, "fasta"):
			headers.append(record.description)
	
	diamondf = pd.read_csv(diamondout, index_col=0)
	viphogdf = pd.read_csv(viphogsout, delimiter="\t", index_col=0)
	phagcndf = pd.read_csv(diamondout, index_col=0)
	genomdf = pd.read_csv(genomadout, index_col=0)
	
	diamondf = diamondf.add_prefix("diamond_")
	viphogdf = viphogdf.add_prefix("viphogs_")
	phagcndf = phagcndf.add_prefix("phagcn_")
	genomdf = genomdf.add_prefix("genomad_")
	
	
	mergedf = pd.DataFrame(index=headers)
	
	mergedf = mergedf.merge(diamondf, left_index=True, right_index=True, how='outer')
	mergedf = mergedf.merge(viphogdf, left_index=True, right_index=True, how='outer')
	mergedf = mergedf.merge(phagcndf, left_index=True, right_index=True, how='outer')
	mergedf = mergedf.merge(genomdf, left_index=True, right_index=True, how='outer')
	
	mergedf = mergedf.replace(["Na", "Unclassified"], np.nan)
	mergedf = mergedf.drop(columns=mergedf.columns[mergedf.columns.str.contains('subfamily', case=False)])
	
	
	mergedf.to_csv(out_csv)
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.description = "Takes the output of different taxonomic annotation tools and returns a merged dataframe in csv format with column prefixes."
	parser.add_argument("--diamondout", required=True, type=str, help="Diamond final taxonomic annotation csv file")
	parser.add_argument("--viphogsout", required=True, type=str, help="Viphogs final taxonomic annotation csv fille")
	parser.add_argument("--phagcnout", required=True, type=str, help="PhaGCN final taxonomic annotation csv file")
	parser.add_argument("--genomadout", required=True, type=str, help="geNomad final taxonomic annotation csv file")
	parser.add_argument("--contigs", required=True, type=str, help="Contigs fasta file used as input for taxonomic annotation")
	parser.add_argument("--output", required=True, type=str, help="Path to merged output file csv")

	args = parser.parse_args()

	merge_taxonomy(args.diamondout, args.viphogsout, args.phagcnout, args.genomadout, args,contigs, args.output)

