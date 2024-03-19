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

echo "Extracting spacers and 3' extensions from R1 fastq files:"
for filename in *R1_001.fastq.gz; do
	filename_array=(${filename//_/ })
	shortname=${filename_array[0]}"_R1_001"
	cutadapt -j 0 -g GAAAAGCCTTGTTTG --discard-untrimmed -o $output_dir$shortname"_5trim_temp.fastq.gz" $filename >> $logfile
	cutadapt -j 0 -a GTTTAAGAGCTATGC --no-indels -m 20 -M 20 --discard-untrimmed -o $output_dir$shortname"_protospacer.fastq.gz" $output_dir$shortname"_5trim_temp.fastq.gz" >> $logfile
	cutadapt -j 0 -g GCACCGAGTCGGTGC --no-indels --discard-untrimmed -o $output_dir$shortname"_5trim_temp.fastq.gz" $filename >> $logfile
	cutadapt -j 0 -a CGCGGTTCTATCTAG --no-indels --discard-untrimmed -o $output_dir$shortname"_RTT-PBS.fastq.gz" $output_dir$shortname"_5trim_temp.fastq.gz" >> $logfile
	echo "$filename processed"
done

echo "Extracting surrogate targets and BCs from R2 fastq files:"
for filename in *R2_001.fastq.gz; do
	filename_array=(${filename//_/ })
	shortname=${filename_array[0]}"_R2_001"
	seqkit seq -r -p -v -t "dna" $filename | gzip -c > $output_dir$shortname"_rc_temp.fastq.gz"
	cutadapt -j 0 -a GGTAAGCTCGTACCG --discard-untrimmed -o $output_dir$shortname"_sensor-BCs_temp.fastq.gz" $output_dir$shortname"_rc_temp.fastq.gz" >> $logfile
	cutadapt -j 0 -l -10 -o $output_dir$shortname"_library-BC.fastq.gz" $output_dir$shortname"_sensor-BCs_temp.fastq.gz" >> $logfile
	cutadapt -j 0 -l -26 -o $output_dir$shortname"_pegBC_temp.fastq.gz" $output_dir$shortname"_sensor-BCs_temp.fastq.gz" >> $logfile
	cutadapt -j 0 -l 16 -o $output_dir$shortname"_pegRNA-BC.fastq.gz" $output_dir$shortname"_pegBC_temp.fastq.gz" >> $logfile
	cutadapt -j 0 -l -81 -o $output_dir$shortname"_sensor_temp.fastq.gz" $output_dir$shortname"_sensor-BCs_temp.fastq.gz" >> $logfile
	cutadapt -j 0 -l 55 -o $output_dir$shortname"_sensor.fastq.gz" $output_dir$shortname"_sensor_temp.fastq.gz" >> $logfile
	echo "$filename processed"
done
rm -- "$output_dir"*_temp.fastq.gz

echo "pegRNA features successfully extracted from fastq files and written to: $script_dir/$output_dir"
cd $work_dir

exit 0
