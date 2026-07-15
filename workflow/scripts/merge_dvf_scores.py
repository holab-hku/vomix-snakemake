import sys
import pandas as pd


def main():
    # The last argument is the output file
    output_file = sys.argv[-1]
    # All preceding arguments are the input files
    input_files = sys.argv[1:-1]

    # Read all files into a list of dataframes
    dfs = [pd.read_csv(f, sep="\t") for f in input_files]

    # Concatenate them safely
    merged_df = pd.concat(dfs, ignore_index=True)

    # Write to output
    merged_df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    main()
