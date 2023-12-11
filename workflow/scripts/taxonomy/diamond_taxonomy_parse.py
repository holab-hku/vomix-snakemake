import pandas as pd 


# name of taxonomy columns in ncbi tax df 
# will be used later for tax assignments
taxcols = ['Species', 'Genus', 'Family', 'Molecule_type', 'Organism_Name', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']

dmndtsvout="results/taxonomy/viral/intermediate/diamond/diamond_out.tsv"
ncbitaxcsv="workflow/database/ncbi/ncbi-virus/output.csv"

taxdf = pd.read_csv(ncbitaxcsv, delimiter=",")
taxdf['Accession_nv'] = taxdf['Accession'].str.rsplit('.', n=1).apply(lambda x: x[0])

dmndheaders = [
		"query_id", 
		"subject_id", 
		"percentage_match", 
		"ali_length", 
		"n_mistmatch", 
		"n_gaps", 
		"start_query", 
		"end_query", 
		"start_subject", 
		"end_subject", 
		"eval", 
		"bit_score" 
		]
df = pd.read_csv(dmndtsvout, delimiter="\t", comment='#', names=dmndheaders)

# add extra columns
df['database_id'] = df['subject_id'].str.rsplit('.', n=1).apply(lambda x: x[0])
df['contig_id'] = df['query_id'].str.rsplit('_', n=1).apply(lambda x: x[0])

mergedf=df.merge(taxdf, left_on="database_id", right_on="Accession_nv")

# group by query id and sort by bit-score
groupdf = mergedf.sort_values(['query_id', 'bit_score'], ascending=False).groupby('query_id')

# 1) FILTER TOP HITS PER PROTEIN
# per protein diamond hit, pick the highest score hit 
# if it does not have taxonomy annotation, pick the second if it is within 
# 25% of the top hit and has taxonomy annotation 

filterlist = []
for group_n, group_df in groupdf:
	# check if highest bit-score has taxonomy
	firstrow = group_df.iloc[0]
	tax = firstrow.loc[taxcols]
	notax = tax.isnull().all()
	nhits = group_df.shape[0]
	# if there is taxonomy or only one top hit
	if (not notax) or (nhits == 1):
		filterlist.append(firstrow)
	else:
		secondrow = group_df.iloc[1]
		tax = secondrow.loc[taxcols]
		notax = tax.isnull().all()
		firsthitbit = firstrow['bit_score']
		secondhitbit = secondrow['bit_score']
		if (secondhitbit * 1.25 >= firsthitbit) and (not notax):
			filterlist.append(secondrow)
		else:
			filterlist.append(firstrow)

tophitsdf = pd.DataFrame(filterlist)

# 2) AGGREGATE PER CONTIG & BIT-SCORE PER TAX RANK
# if the sum of the bit score for a species is 70% of the total 
# accumulated bitscore [since bitscores are independent of query length]
# the contig will be assigned that value

taxdict = {}
# Iterate through each contig 
for contig in tophitsdf['contig_id'].unique():
	taxdict[contig] = {}
	# Iterate through each tax level 
	for taxlev in taxcols:
		contigdf = tophitsdf[tophitsdf['contig_id'].values == contig]
		groupdf = contigdf.groupby(taxlev)
		bitsumdf = groupdf['bit_score'].sum()
		bitsumcontig = bitsumdf.sum()
		bitpercdf = pd.DataFrame(bitsumdf / bitsumcontig * 100).reset_index()
		bitpercdf  = bitpercdf.rename(columns={'bit_score':'bit_perc'})
		# if the function fails it mean there is not tax annotation at that level
		try:
			maxhit = bitpercdf.iloc[bitpercdf['bit_perc'].idxmax()]
		except ValueError:
			tax = "Na"
			taxdict[contig][taxlev] = tax
			continue
		score = maxhit['bit_perc']
		if score > 0.7:
			tax = maxhit[taxlev]
			taxdict[contig][taxlev] = tax
		else:
			tax = "Unclassified"
			taxdict[contig][taxlev] = tax
			continue

taxdf = pd.DataFrame(taxdict).T

taxdf.to_csv("testout.csv")
