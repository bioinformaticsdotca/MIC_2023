#!/bin/bash

# Get the conda environment path at the start of the shell script.
source ~/.bashrc

# Activate the seqkit conda environment.
conda activate seqkit_env

# The path to a list of sample ids.
ids_infile="KGHS_pilot_subset_4_sample_list.txt"

# A minimum length of 1500 bps. 
minimum_length=1500

# The output directory to write the filtered files.
output_dir="filtered_metagenomes"

# Create the output directory if it does not exist.
mkdir -p $output_dir

# Internal Field Separator (IFS) used when using cat in a for loop so that lines are separated by newlines instead of separated based on spaces.
IFS=$'\n'

# Iterate over the list of sample ids 
for sample_id in $(cat $ids_infile); 
do
	# The sample name output directory.
	sample_dir="${output_dir}/${sample_id}"

    # Create the sample name output directory if it does not exist.
    mkdir -p $sample_dir

	# Run the seqkit command.
	echo "seqkit seq -m ${minimum_length} metagenome_assemblies/${sample_id}/scaffolds.fasta > ${sample_dir}/${sample_id}_min${minimum_length}.fasta"
	seqkit seq -m ${minimum_length} metagenome_assemblies/${sample_id}/scaffolds.fasta > ${sample_dir}/${sample_id}_min${minimum_length}.fasta

done




