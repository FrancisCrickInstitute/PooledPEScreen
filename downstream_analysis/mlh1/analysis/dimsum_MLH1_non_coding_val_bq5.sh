# Prime editing screen S002, validation of hits in 4 non coding regions of the MLH1 gene, sequencing of the endogenous site, (2x upstream, x11, i15): Run dimsum pipeline with per base quality setting = 5 on:
# all D20/D34 combinations
# all D20/negative ctrl combinations

input_dir="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/240129_Screen002_val/fastq"
output_dir="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/240129_Screen002_val/output/dimsum_output"

#---------------------------------------------------------------------------------------------------------------------------------
# upstream region 1 - D20/D34

first5adapter="GCCAGATCACCTCAGCAGAG"
second5adapter="TGGACAGCTTGAATGCCAGT"
wildtype="GCACACAAGCCCGGTTCCGGCATCTCTGCTCCTATTGGCTGGATATTTCGTATTCCCCGAGCTCCTAAAAACGAACCAATAGGAAGAGCGGACAGCGATCTCTAACGCGCAAGCGCATATCCTTCTAGGTAGCGGGCAGTAGCCGCTTCAGGGAGGGACGAAGAGACCCAGCAACCCACAGAGTTGAGAAATTTG"
max_errors="0.5"
max_substitutions="3"
CPU_cores="4"

base_quality="5"

experiment_design="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/240129_Screen002_val/input/dimsum_experimental_design_sheet/CK_240129_S002_val_dimsum_experimental_design_u1_D20_D34.txt"
output_name="CK_20240129_dimsum_u1_D20_D34_bq5"

DiMSum --fastqFileDir $input_dir --experimentDesignPath $experiment_design --outputPath $output_dir --projectName $output_name --cutadapt5First $first5adapter --cutadapt5Second $second5adapter --wildtypeSequence $wildtype --vsearchMinQual $base_quality --vsearchMaxee $max_errors --maxSubstitutions $max_substitutions --retainIntermediateFiles T --numCores $CPU_cores

# upstream region 1 - neg/D20

first5adapter="GCCAGATCACCTCAGCAGAG"
second5adapter="TGGACAGCTTGAATGCCAGT"
wildtype="GCACACAAGCCCGGTTCCGGCATCTCTGCTCCTATTGGCTGGATATTTCGTATTCCCCGAGCTCCTAAAAACGAACCAATAGGAAGAGCGGACAGCGATCTCTAACGCGCAAGCGCATATCCTTCTAGGTAGCGGGCAGTAGCCGCTTCAGGGAGGGACGAAGAGACCCAGCAACCCACAGAGTTGAGAAATTTG"
max_errors="0.5"
max_substitutions="3"
CPU_cores="4"

base_quality="5"

experiment_design="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/240129_Screen002_val/input/dimsum_experimental_design_sheet/CK_240129_S002_val_dimsum_experimental_design_u1_neg.txt"
output_name="CK_20240129_dimsum_u1_neg_bq5"

DiMSum --fastqFileDir $input_dir --experimentDesignPath $experiment_design --outputPath $output_dir --projectName $output_name --cutadapt5First $first5adapter --cutadapt5Second $second5adapter --wildtypeSequence $wildtype --vsearchMinQual $base_quality --vsearchMaxee $max_errors --maxSubstitutions $max_substitutions --retainIntermediateFiles T --numCores $CPU_cores

#---------------------------------------------------------------------------------------------------------------------------------
# upstream region 2 - D20/D34

first5adapter="AGCTGTCCAATCAATAGCTGC"
second5adapter="CGATGCGGTTCACCACTGT"
wildtype="CGCTGAAGGGTGGGGCTGGATGGCGTAAGCTACAGCTGAAGGAAGAACGTGAGCACGAGGCACTGAGGTGATTGGCTGAAGGCACTTCCGTTGAGCATCTAGACGTTTCCTTGGCTCTTCTGGCGCCAAAATGTCGTTCGTGGCAGGGGTTATTCGGCGGCTGGACGAG"
max_errors="0.5"
max_substitutions="3"
CPU_cores="4"

base_quality="5"

experiment_design="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/240129_Screen002_val/input/dimsum_experimental_design_sheet/CK_240129_S002_val_dimsum_experimental_design_u2_D20_D34.txt"
output_name="CK_20240129_dimsum_u2_D20_D34_bq5"

DiMSum --fastqFileDir $input_dir --experimentDesignPath $experiment_design --outputPath $output_dir --projectName $output_name --cutadapt5First $first5adapter --cutadapt5Second $second5adapter --wildtypeSequence $wildtype --vsearchMinQual $base_quality --vsearchMaxee $max_errors --maxSubstitutions $max_substitutions --retainIntermediateFiles T --numCores $CPU_cores

# upstream region 2 - neg/D20

first5adapter="AGCTGTCCAATCAATAGCTGC"
second5adapter="CGATGCGGTTCACCACTGT"
wildtype="CGCTGAAGGGTGGGGCTGGATGGCGTAAGCTACAGCTGAAGGAAGAACGTGAGCACGAGGCACTGAGGTGATTGGCTGAAGGCACTTCCGTTGAGCATCTAGACGTTTCCTTGGCTCTTCTGGCGCCAAAATGTCGTTCGTGGCAGGGGTTATTCGGCGGCTGGACGAG"
max_errors="0.5"
max_substitutions="3"
CPU_cores="4"

base_quality="5"

experiment_design="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/240129_Screen002_val/input/dimsum_experimental_design_sheet/CK_240129_S002_val_dimsum_experimental_design_u2_neg.txt"
output_name="CK_20240129_dimsum_u2_neg_bq5"

DiMSum --fastqFileDir $input_dir --experimentDesignPath $experiment_design --outputPath $output_dir --projectName $output_name --cutadapt5First $first5adapter --cutadapt5Second $second5adapter --wildtypeSequence $wildtype --vsearchMinQual $base_quality --vsearchMaxee $max_errors --maxSubstitutions $max_substitutions --retainIntermediateFiles T --numCores $CPU_cores

#---------------------------------------------------------------------------------------------------------------------------------
# x11 region - D20/D34

first5adapter="CCCCCTCCCACTATCTAAGGTAAT"
second5adapter="AAGGCCCCAGAGAAGTAGCT"
wildtype="TGTTCTCTCTTATTTTCCTGACAGTTTAGAAATCAGTCCCCAGAATGTGGATGTTAATGTGCACCCCACAAAGCATGAAGTTCACTTCCTGCACGAGGAGAGCATCCTGGAGCGGGTGCAGCAGCACATCGAGAGCAAGCTCCTGGGCTCCAATTCCTCCAGGATGTACTTCACCCAGGTCAGGGCGCTTCTCATCC"
max_errors="0.5"
max_substitutions="3"
CPU_cores="4"

base_quality="5"

experiment_design="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/240129_Screen002_val/input/dimsum_experimental_design_sheet/CK_240129_S002_val_dimsum_experimental_design_x11_D20_D34.txt"
output_name="CK_20240129_dimsum_x11_D20_D34_bq5"

DiMSum --fastqFileDir $input_dir --experimentDesignPath $experiment_design --outputPath $output_dir --projectName $output_name --cutadapt5First $first5adapter --cutadapt5Second $second5adapter --wildtypeSequence $wildtype --vsearchMinQual $base_quality --vsearchMaxee $max_errors --maxSubstitutions $max_substitutions --retainIntermediateFiles T --numCores $CPU_cores

# x11 region - neg/D20

first5adapter="CCCCCTCCCACTATCTAAGGTAAT"
second5adapter="AAGGCCCCAGAGAAGTAGCT"
wildtype="TGTTCTCTCTTATTTTCCTGACAGTTTAGAAATCAGTCCCCAGAATGTGGATGTTAATGTGCACCCCACAAAGCATGAAGTTCACTTCCTGCACGAGGAGAGCATCCTGGAGCGGGTGCAGCAGCACATCGAGAGCAAGCTCCTGGGCTCCAATTCCTCCAGGATGTACTTCACCCAGGTCAGGGCGCTTCTCATCC"
max_errors="0.5"
max_substitutions="3"
CPU_cores="4"

base_quality="5"

experiment_design="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/240129_Screen002_val/input/dimsum_experimental_design_sheet/CK_240129_S002_val_dimsum_experimental_design_x11_neg.txt"
output_name="CK_20240129_dimsum_x11_neg_bq5"

DiMSum --fastqFileDir $input_dir --experimentDesignPath $experiment_design --outputPath $output_dir --projectName $output_name --cutadapt5First $first5adapter --cutadapt5Second $second5adapter --wildtypeSequence $wildtype --vsearchMinQual $base_quality --vsearchMaxee $max_errors --maxSubstitutions $max_substitutions --retainIntermediateFiles T --numCores $CPU_cores

#---------------------------------------------------------------------------------------------------------------------------------
# i15 region - D20/D34

first5adapter="TGATCTGCCTGTGGTTTTCT"
second5adapter="ACTCATCCCCCATATACCCCA"
wildtype="AGAAATACGAAAGTTGAGTCCTTAAGGCTACACAGAAAGAAAGTACCTCCCCAGGGCTTCACCCTTCCCATCCTTTCAGCAGGCTTTTTGTCTGTCGTATCTTCTCTGTTGAAATGGCCATTGACAAGAGGAGGAAAGGGGTTTTGTTGTGGATTGTTCAGGCACTTCCTT"
max_errors="0.5"
max_substitutions="3"
CPU_cores="4"

base_quality="5"

experiment_design="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/240129_Screen002_val/input/dimsum_experimental_design_sheet/CK_240129_S002_val_dimsum_experimental_design_i15_D20_D34.txt"
output_name="CK_20240129_dimsum_i15_D20_D34_bq5"

DiMSum --fastqFileDir $input_dir --experimentDesignPath $experiment_design --outputPath $output_dir --projectName $output_name --cutadapt5First $first5adapter --cutadapt5Second $second5adapter --wildtypeSequence $wildtype --vsearchMinQual $base_quality --vsearchMaxee $max_errors --maxSubstitutions $max_substitutions --retainIntermediateFiles T --numCores $CPU_cores

# i15 region - neg/D20

first5adapter="TGATCTGCCTGTGGTTTTCT"
second5adapter="ACTCATCCCCCATATACCCCA"
wildtype="AGAAATACGAAAGTTGAGTCCTTAAGGCTACACAGAAAGAAAGTACCTCCCCAGGGCTTCACCCTTCCCATCCTTTCAGCAGGCTTTTTGTCTGTCGTATCTTCTCTGTTGAAATGGCCATTGACAAGAGGAGGAAAGGGGTTTTGTTGTGGATTGTTCAGGCACTTCCTT"
max_errors="0.5"
max_substitutions="3"
CPU_cores="4"

base_quality="5"

experiment_design="/camp/home/kajbac/home/shared/projects/PE/ScreeningData/240129_Screen002_val/input/dimsum_experimental_design_sheet/CK_240129_S002_val_dimsum_experimental_design_i15_neg.txt"
output_name="CK_20240129_dimsum_i15_neg_bq5"

DiMSum --fastqFileDir $input_dir --experimentDesignPath $experiment_design --outputPath $output_dir --projectName $output_name --cutadapt5First $first5adapter --cutadapt5Second $second5adapter --wildtypeSequence $wildtype --vsearchMinQual $base_quality --vsearchMaxee $max_errors --maxSubstitutions $max_substitutions --retainIntermediateFiles T --numCores $CPU_cores