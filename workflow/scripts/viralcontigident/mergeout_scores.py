import sys
import os
import pandas as pd


def merge_csv_files(files):
	merged_data = pd.DataFrame()
	data_frames = []
	for f in files:
		df = pd.read_csv(f, index_col=0)
		sample_id = f.split("/")[3]
		df['sample_id'] = sample_id  # Add a new column with the file name
		df.index = df.index.astype(str) + "_" + str(sample_id)
		data_frames.append(df)

	merged_data = pd.concat(data_frames, ignore_index=False)
	# Print merged table to STDOUT
	print(merged_data.to_csv(index=True))

# Retrieve the list of files from command-line arguments
csv_files = sys.argv[1:]

# Check if at least one file is provided
if len(csv_files) < 1:
	print("Please provide at least one CSV file.")
else:
	merge_csv_files(csv_files)

