"""
    Project :  Prime editing pilot screens. Revision experiment. Script to convert dimsum output RD data to .csv files and save the "all_variants" dataframe.
    Date : 240909
    Python version 3.10

"""

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import pandas as pd
import os
import sys

# Activate automatic conversion between R and Pandas
pandas2ri.activate()


def rdata_to_csv(rdata_file, output_dir, R_df_name):
    # Load the RData file
    ro.r["load"](rdata_file)

    # Check if the desired data frame exists in the R environment
    if R_df_name in ro.r.ls():
        # Get the desired data frame object
        r_dataframe = ro.r[R_df_name]

        # Check if the object is a data frame
        if "data.frame" in ro.r["class"](r_dataframe):
            # Convert the R data frame to a Pandas DataFrame
            df = pandas2ri.rpy2py(r_dataframe)

            # Extract the base name of the .RData file and create the .csv filename
            base_name = os.path.splitext(os.path.basename(rdata_file))[0]
            output_csv = os.path.join(output_dir, f"{base_name}.csv")

            # Save the DataFrame as a CSV file
            df.to_csv(output_csv, index=False)
            print(f"Saved {R_df_name} from {rdata_file} to {output_csv}")
        else:
            print(f"Skipping {R_df_name} in {rdata_file}, not a data frame.")
    else:
        print(f"{R_df_name} not found in {rdata_file}.")


def process_folder(input_folder, output_folder, R_df_name):
    # Check if output directory exists, if not, create it
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Loop through all .RData files in the input folder
    for filename in os.listdir(input_folder):
        if filename.endswith(".RData"):
            rdata_file_path = os.path.join(input_folder, filename)
            rdata_to_csv(rdata_file_path, output_folder, R_df_name)


if __name__ == "__main__":
    # Check if correct number of arguments are provided
    if len(sys.argv) != 4:
        print(
            "Usage: python rdata_to_csv_specific_df.py <input_folder> <output_folder> <specific_dataframe_name>"
        )
        sys.exit(1)

    # Get the input folder, output folder, and the specific data frame name
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    R_df_name = sys.argv[3]

    # Process all .RData files in the input folder
    process_folder(input_folder, output_folder, R_df_name)
