# Prime editing screen S002 MLH1 endogenous site: Run dimsum pipeline with per base quality setting = 5 
# on all negative control/D20 combinations
# include ROI adjacent region for identifying errors



input_dir="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/230807_Screen002/fastqs/MLH1/endogenous/"
output_dir="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/230807_Screen002/output/MLH1/dimsum_output"

first5adapter="TTCATGACTTTGTGTGAATGTACA"
second5adapter="AGCATGCTCATCTCTTTCAAAG"
wildtype="CCTGTGACCTCACCCCTCAGGACAGTTTTGAACTGGTTGCTTTCTTTTTATTGTTTAGATCGTCTGGTAGAATCAACTTCCTTGAGAAAAGCCATAGAAACAGTGTATGCAGCCTATTTGCCCAAAAACACACACCCATTCCTGTACCTCAGGTAATGTAGCACCAAACTCCTCAACCAAGACTCACAAGGAACAGATGTTCTATCAGGCTCTCCT"
max_errors="0.5"
max_substitutions="3"
CPU_cores="4"

base_quality="5"

#----------------------------------------------------------------------------------------------------------------------------------------------
# Run dimsum on all samples with per base quality set to 5, neg/D20

# S1 negative control compared to D20
experiment_design="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/230807_Screen002/input/MLH1/dimsum_experimental_design_sheet/CK_231002_S002_dimsum_experimental_design_S1_neg_D20.txt"
output_name="CK_20231002_dimsum_S1_neg_D20_bq5"

DiMSum --fastqFileDir $input_dir --experimentDesignPath $experiment_design --outputPath $output_dir --projectName $output_name --cutadapt5First $first5adapter --cutadapt5Second $second5adapter --wildtypeSequence $wildtype --vsearchMinQual $base_quality --vsearchMaxee $max_errors --maxSubstitutions $max_substitutions --retainIntermediateFiles T --numCores $CPU_cores

# S2 negative control compared to D20
experiment_design="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/230807_Screen002/input/MLH1/dimsum_experimental_design_sheet/CK_231002_S002_dimsum_experimental_design_S2_neg_D20.txt"
output_name="CK_20231002_dimsum_S2_neg_D20_bq5"

DiMSum --fastqFileDir $input_dir --experimentDesignPath $experiment_design --outputPath $output_dir --projectName $output_name --cutadapt5First $first5adapter --cutadapt5Second $second5adapter --wildtypeSequence $wildtype --vsearchMinQual $base_quality --vsearchMaxee $max_errors --maxSubstitutions $max_substitutions --retainIntermediateFiles T --numCores $CPU_cores

# S1 ouabain selected negative control compared to D20
experiment_design="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/230807_Screen002/input/MLH1/dimsum_experimental_design_sheet/CK_231002_S002_dimsum_experimental_design_S1_O_neg_D20.txt"
output_name="CK_20231002_dimsum_S1_O_neg_D20_bq5"

DiMSum --fastqFileDir $input_dir --experimentDesignPath $experiment_design --outputPath $output_dir --projectName $output_name --cutadapt5First $first5adapter --cutadapt5Second $second5adapter --wildtypeSequence $wildtype --vsearchMinQual $base_quality --vsearchMaxee $max_errors --maxSubstitutions $max_substitutions --retainIntermediateFiles T --numCores $CPU_cores

# S2 ouabain selected negative control compared to D20
experiment_design="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/230807_Screen002/input/MLH1/dimsum_experimental_design_sheet/CK_231002_S002_dimsum_experimental_design_S2_O_neg_D20.txt"
output_name="CK_20231002_dimsum_S2_O_neg_D20_bq5"

DiMSum --fastqFileDir $input_dir --experimentDesignPath $experiment_design --outputPath $output_dir --projectName $output_name --cutadapt5First $first5adapter --cutadapt5Second $second5adapter --wildtypeSequence $wildtype --vsearchMinQual $base_quality --vsearchMaxee $max_errors --maxSubstitutions $max_substitutions --retainIntermediateFiles T --numCores $CPU_cores