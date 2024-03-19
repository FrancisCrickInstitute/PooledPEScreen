#!/bin/bash

# Execute from within directory containing input fastq files by calling:
# (chmod +x extract_features_from_fastqs.sh)
# time ./extract_features_from_fastqs.sh

work_dir=$(pwd)
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
cd $script_dir
output_dir=processed/
mkdir -p $output_dir
datum=$(date +"%Y%m%d")
logfile=$output_dir$datum"_cutadaptlog.txt"
touch $logfile

echo "Extracting spacer and scaffold from R1 fastq files:"
for filename in *R1_001.fastq.gz; do
	filename_array=(${filename//_/ })
	shortname=${filename_array[0]}"_R1_001"
	cutadapt -j 0 -a "^ATGCGAGAAAAGCCTTGTTTG...GTTTHAGAGCTA" --discard-untrimmed -o $output_dir$shortname"_spacer.fastq.gz" $filename >> $logfile
	cutadapt -j 0 -a "^ATGCGAGAAAAGCCTTGTTTGN{20}GTTT...GCACCGAGTCGGTGC" --discard-untrimmed -o $output_dir$shortname"_scaffold.fastq.gz" $filename >> $logfile
	echo "$filename processed"
done

echo "Spacer and scaffold successfully extracted from fastq files and written to: $script_dir/$output_dir"
cd $work_dir

exit 0
