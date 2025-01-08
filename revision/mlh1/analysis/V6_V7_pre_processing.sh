#Date: 04/11/24
#Experiment: PE manuscript revision experiment
#Received de-multiplexed fastq files from ASF.
#Run Seqprep pipeline for QC and read merging. 

# seq_dir="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/240822_Revision/validation_of_splice_variants/fastqs/MLH1"
out_dir="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/240822_Revision/output/"
adapter1="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC" #this is the nextera reverse adapter, as seen in R1
adapter2="CTGTCTCTTATACACATCTGACGCTGCCGACGA"  #this is the nextera forward adapter, as seen in R2
seq_type="M" #"N" for nextseq, "M" for Miseq

ml load Python/2.7.18-GCCcore-9.3.0

#count how many times "@M" is found in each fastq, will be once per read for a Miseq (change seqtype if Nextseq run)
zgrep -c @$seq_type *R1_* | tee read_counts_R1.txt

zgrep -c @$seq_type *R2_* | tee read_counts_R2.txt


mkdir Seqprep
mkdir Seqprep/R1
mkdir Seqprep/R2
mkdir Seqprep/merged
python ~/home/users/findlag/bin/run_seqprep_VHL_pipeline.py $adapter1 $adapter2
echo "Running Seqprep on all samples in parallel."
sh run_seqprep.sh
cd Seqprep/merged
echo "Seqprep done."
zgrep -c @$seq_type *.fastq.gz >> seqprep_read_counts.txt

mkdir no_Ns
python ~/home/users/findlag/bin/run_remove_n_bases.py #operates with getcwd() uses remove_n_bases.py script - this step seems only necessary for some NS runs
sh run_remove_n_bases.sh
echo "remove_n_bases done."
cd no_Ns
zgrep -c @$seq_type *.fastq >> no_Ns_read_counts.txt

