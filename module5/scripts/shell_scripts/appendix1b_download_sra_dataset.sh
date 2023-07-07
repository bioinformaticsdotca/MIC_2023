#!/bin/bash

# Get the conda environment path at the start of the shell script.
source ~/.bashrc

# Activate the NCBI SRA tools  conda environment.
conda activate sra_tools_env

# The path to a list of the sample ids we want to grab from the SRA Run Info file.
ids_infile="KGHS_pilot_subset_4_sample_list.txt"

# The SRA run information summary file.
sra_run_info_infile="SraRunInfo.csv"

# The output directory to download SRA files.
fastq_output_dir="cleaned_fastq_files"

# The temporary file directory for fasterq-dump.
tmp_dir="$fastq_output_dir/tmp_dir"

# The number of threads to use for fasterq-dump.
num_threads=8

# Internal Field Separator (IFS) used when using cat in a for loop so that lines are separated by newlines instead of separated based on spaces.
IFS=$'\n'

# Iterate over the SraRunInfo.csv file 
for sra_run_info_entry in $(tail -n+2 $sra_run_info_infile); 
do 
	echo $sra_run_info_entry;

	# Get the SRA Run ID from the sra_run_info_infile
	sra_run_id=$(echo $sra_run_info_entry | cut -d ',' -f1); 
	echo $sra_run_id;

	# Get the sample_id from the sra_run_info_infile row.
	sample_id=$(echo $sra_run_info_entry | cut -d ',' -f30)
		
	# Download the SRA fastq files for each sample_id using the sra_run_id SRA Run ID. Split them into read1 and read2 fastqs. Use num_threads for the number of threads for fasterq-dump. 
	echo "fasterq-dump --threads $num_threads --outfile $sample_id --outdir ${fastq_output_dir} --temp $tmp_dir --split-files $sra_run_id"
	fasterq-dump --threads $num_threads --outfile $sample_id --outdir ${fastq_output_dir} --temp $tmp_dir --split-files $sra_run_id

done




