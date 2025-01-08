"""
    Project :  Prime editing pilot screens. Revision experiment. Script tp generate input txt files (experimental design sheets) for the dimsum pipeline. 
    Date : 240906
    Python version 3.10

"""
# IMPORT
import os
import pandas as pd
import pathlib as pl

# PATHS

fastq_directory = "/camp/home/kajbac/home/users/kajbac/fastq/240822"

# CSV file specifying sample comparrisons
dimsum_map_path = "/camp/home/kajbac/home/shared/projects/PE/ScreeningData/240822_Revision/input/dimsum_map/CK_dimsum_map.csv"
dimsum_map = pd.read_csv(dimsum_map_path)
output_path = pl.Path(
    "/camp/home/kajbac/home/shared/projects/PE/ScreeningData/240822_Revision/input/dimsum_experimental_design_sheet"
)


# match R1 and R2 filename pairs for a given identifier
def find_file_pairs(identifier, directory):
    pair1 = pair2 = None
    for file_name in os.listdir(directory):
        if identifier in file_name:
            if "R1" in file_name:
                pair1 = file_name
            elif "R2" in file_name:
                pair2 = file_name
    return pair1, pair2


# Iterate over each row of the dimsum_map
for index, row in dimsum_map.iterrows():
    sample_name = row["Sample"]
    input_str = row["input"]
    output_str = row["output"]

    # Create empty experimental design txt file
    output_file = output_path / f"CK_{sample_name}.txt"

    # Find file pairs for D19 and D33
    input_pair1, input_pair2 = find_file_pairs(input_str, fastq_directory)
    output_pair1, output_pair2 = find_file_pairs(output_str, fastq_directory)

    # Open the experimental design sheet to store the file name pairs
    with open(output_file, "w") as f:
        # Write the header
        f.write(
            "sample_name\texperiment\tselection_id\tselection_replicate\ttechnical_replicate\tpair1\tpair2\n"
        )

        # Write the 'input' row
        f.write(f"input\t1\t0\t\t\t{input_pair1}\t{input_pair2}\n")

        # Write the 'output' row
        f.write(f"output\t1\t1\t1\t\t{output_pair1}\t{output_pair2}\n")

    print(f"Created {sample_name} with paired files for input and output.")
