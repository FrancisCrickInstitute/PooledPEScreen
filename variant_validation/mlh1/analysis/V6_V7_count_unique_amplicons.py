"""
    Project :  Prime editing pilot screens. Revisions. V6 and V7 splicing outcome analysis. Script to count unique sequences in fastq files and generate sequence count output files. 
    Date : 241115
    Python version 3.10
"""
# IMPORT

import os
import argparse
from collections import Counter
import pandas as pd


# FUNCTIONS
# Process a single FASTQ file
def process_fastq(file_path):
    sequences = []

    with open(file_path, "r") as f:
        for i, line in enumerate(f):
            # Sequences are typically every 4th line in FASTQ files
            if i % 4 == 1:
                sequences.append(line.strip())

    # Count occurrences of each unique sequence
    sequence_counts = Counter(sequences)

    # Create dataframe
    df = pd.DataFrame(
        {
            "Sequence": list(sequence_counts.keys()),
            "Count": list(sequence_counts.values()),
        }
    )

    return df


def main():
    # Argument parser for command line inputs
    parser = argparse.ArgumentParser(
        description="Process FASTQ files and generate sequence counts."
    )
    parser.add_argument(
        "input_path", help="Path to the FASTQ file or directory containing FASTQ files."
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        help="Directory to save output CSV files. Defaults to the same directory as the input.",
        default=None,
    )
    args = parser.parse_args()

    input_path = args.input_path
    output_dir = args.output_dir if args.output_dir else os.path.dirname(input_path)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if os.path.isfile(input_path):
        # Process a single file
        file_name = os.path.basename(input_path)
        df = process_fastq(input_path)
        output_file = os.path.join(
            output_dir, file_name.replace(".fastq", "_counts.csv")
        )
        df.to_csv(output_file, index=False)
        print(f"Processed {file_name} and saved to {output_file}")

    elif os.path.isdir(input_path):
        # Process all FASTQ files in a directory
        for file_name in os.listdir(input_path):
            if file_name.endswith(".fastq"):
                file_path = os.path.join(input_path, file_name)
                print(f"Processing {file_name}...")
                df = process_fastq(file_path)
                output_file = os.path.join(
                    output_dir, file_name.replace(".fastq", "_counts.csv")
                )
                df.to_csv(output_file, index=False)
                print(f"Saved dataframe to {output_file}.")
    else:
        print(f"Error: {input_path} is not a valid file or directory.")


if __name__ == "__main__":
    main()
