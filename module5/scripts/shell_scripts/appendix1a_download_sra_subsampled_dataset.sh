#!/bin/bash
#SBATCH --partition=synergy,cpu2013,cpu2019,cpu2021,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=7-00:00:00
#SBATCH --mem=5G
#SBATCH --output=run_download_sra_subsampled_dataset.%A.out
#SBATCH --error=run_download_sra_subsampled_dataset.%A.err

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

# Iterate over the list of SRA ids that we got from the 
for sample_id in $(cat $ids_infile); 
do 
	echo $sample_id;

	# Get the SRA Run ID from the sra_run_info_infile using the sample_id to return the line containing the sample id..	
	sra_run_id=$(grep "$sample_id" < $sra_run_info_infile | cut -d ',' -f1); 
	echo $sra_run_id;
	
	# Download the SRA fastq files for each sample_id using the sra_run_id SRA Run ID. Split them into read1 and read2 fastqs. Use num_threads for the number of threads for fasterq-dump. 
	echo "fasterq-dump --threads $num_threads --outfile $sample_id --outdir ${fastq_output_dir} --temp $tmp_dir --split-files $sra_run_id"
	fasterq-dump --threads $num_threads --outfile $sample_id --outdir ${fastq_output_dir} --temp $tmp_dir --split-files $sra_run_id

done




