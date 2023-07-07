#!/bin/bash

# Get the bashrc information for conda.
source ~/.bashrc

# Activate the spades conda environment to run metaspades.
conda activate spades_env

# The number of cpu threads for metaspades.
num_threads=28

# The RAM memory in GB for metaspades.
memory_in_gb=60

# The list of sample ids.
list_file="KGHS_pilot_subset_4_sample_list.txt"

# The path fastq input file directory.
input_dir="cleaned_fastq_files"

# Internal Field Separator (IFS) used when using cat in a job array so that lines are separated by newlines instead of separated based on spaces.
IFS=$'\n'

for sample_id in $(cat $list_file);
do

	# The read 1 fastq file for each sample.
	fastq_read1="${input_dir}/${sample_id}_1.fastq"

	# The read 2 fastq file for each sample.
	fastq_read2="${input_dir}/${sample_id}_2.fastq"

	# The path metaspades metagenome assembly sample output directory.
	output_dir="metagenome_assemblies/${sample_id}"
	mkdir -p $output_dir

	# Run the metaspades.py program using the read1 and read2 fastq files 
	echo "metaspades.py -k 21,33,55,77,99,127 -t $num_threads -m $memory_in_gb -o $output_dir -1 $fastq_read1 -2 $fastq_read2"
	metaspades.py -k 21,33,55,77,99,127 -t $num_threads -m $memory_in_gb -o $output_dir -1 $fastq_read1 -2 $fastq_read2
done

