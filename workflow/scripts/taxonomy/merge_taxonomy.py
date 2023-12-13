import pandas as pd
from Bio import SeqIO
import numpy as np


#inputs
contigdf

diamondout = "results/taxonomy/viral/intermediate/diamond/taxonomy.csv"
viphogsout = "results/taxonomy/viral/intermediate/viphogs/taxonomy.tsv"
phagcnout = "results/taxonomy/viral/intermediate/phagcn/taxonomy.csv"
genomadout = "results/taxonomy/viral/intermediate/genomad/taxonomy.csv"

contigsfa = "results/viralcontigident/output/combined.final.vOTUs.fa"

#output
outcsv = "test.csv"


#function 
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



taxcols = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']

scoredict = {}
scoredict = {key: {} for key in headers}

for key in scoredict:
	for taxlev in taxcols:
		scoredict[taxlev] = set(mergedf.filter(like=taxlev).values.flatten())





mergedf.to_csv(outcsv)


