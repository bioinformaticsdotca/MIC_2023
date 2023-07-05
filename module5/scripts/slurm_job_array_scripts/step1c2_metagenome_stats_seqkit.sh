#!/bin/bash

# Get the conda environment path at the start of the shell script.
source ~/.bashrc

# Activate the seqkit conda environment.
conda activate seqkit_env

# The path to a list of sample ids.
ids_infile="KGHS_pilot_subset_4_sample_list.txt"

# Make a list of minimum lengths  
minimum_lengths=(1500 2000 2500)

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

	;
	# The assembly scaffold length directory for each sample.
	sample_assembly_len_dir="${sample_dir}/${sample_id}_min${minimum_length}"
		
	# Create the assembly scaffold length directory for each sample.
	mkdir -p $sample_assembly_len_dir

		# Run the seqkit command.
		echo "seqkit stats -a ${sample_assembly_len_dir}/${sample_id}_min${minimum_length}.fasta > ${sample_assembly_len_dir}/${sample_id}_min${minimum_length}.seqkit.stats.txt"
		seqkit stats -a ${sample_assembly_len_dir}/${sample_id}_min${minimum_length}.fasta > ${sample_assembly_len_dir}/${sample_id}_min${minimum_length}.seqkit.stats.txt
	done
done




