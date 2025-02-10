"""
    Project :  Prime editing pilot screens. Revision experiment. Script to generate bash scripts for u2 samples run through the dimsum pipeline. 
    Date : 240906
    Python version 3.10

"""

import os
import argparse
import pathlib as pl


def generate_bash_script(output_file, input_directory, filename):
    # Variable parameters
    base_name = filename.replace(".txt", "")
    experiment_design = f"{input_directory}/{filename}"
    output_name = f"{base_name}_bq5"

    # Fixed parameters for u2
    input_dir = "/camp/home/kajbac/home/users/kajbac/fastq/240822"
    output_dir = "/camp/home/kajbac/home/shared/projects/PE/ScreeningData/240822_Revision/output/dimsum_output"
    first5adapter = "AGCTGTCCAATCAATAGCTGC"
    second5adapter = "CGATGCGGTTCACCACTGT"
    wildtype = "CGCTGAAGGGTGGGGCTGGATGGCGTAAGCTACAGCTGAAGGAAGAACGTGAGCACGAGGCACTGAGGTGATTGGCTGAAGGCACTTCCGTTGAGCATCTAGACGTTTCCTTGGCTCTTCTGGCGCCAAAATGTCGTTCGTGGCAGGGGTTATTCGGCGGCTGGACGAG"
    max_errors = "0.5"
    max_substitutions = "3"
    CPU_cores = "4"

    base_quality = "5"

    # Create Bash script content
    bash_script_content = f"""#!/bin/bash
# This script processes {filename}

input_dir="{input_dir}"
output_dir="{output_dir}"
first5adapter="{first5adapter}"
second5adapter="{second5adapter}"
wildtype="{wildtype}"
max_errors="{max_errors}"
max_substitutions="{max_substitutions}"
CPU_cores="{CPU_cores}"
base_quality="{base_quality}"

experiment_design="{experiment_design}"
output_name="{output_name}"

echo "Experiment design: ${{experiment_design}}"
echo "Output name: ${{output_name}}"

# Run DiMSum with specified parameters
DiMSum --fastqFileDir $input_dir --experimentDesignPath $experiment_design --outputPath $output_dir --projectName $output_name --cutadapt5First $first5adapter --cutadapt5Second $second5adapter --wildtypeSequence $wildtype --vsearchMinQual $base_quality --vsearchMaxee $max_errors --maxSubstitutions $max_substitutions --retainIntermediateFiles T --numCores $CPU_cores
"""

    # Write Bash script content to output file
    with open(output_file, "w") as file:
        file.write(bash_script_content)


def main():
    parser = argparse.ArgumentParser(
        description='Generate Bash scripts for files containing "u2" in their names.'
    )
    parser.add_argument(
        "input_directory",
        type=str,
        help="Directory containing dimsum experimental design sheets",
    )
    parser.add_argument(
        "output_directory",
        type=str,
        help="Directory to save the generated Bash scripts",
    )

    args = parser.parse_args()
    input_directory = args.input_directory
    output_directory = args.output_directory

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Loop through files in the specified directory
    for filename in os.listdir(input_directory):
        if "u2" in filename:
            output_file_path = os.path.join(output_directory, f"{filename}_dimsum.sh")
            generate_bash_script(output_file_path, input_directory, filename)
            print(f"Generated {filename}_dimsum.sh")


if __name__ == "__main__":
    main()
