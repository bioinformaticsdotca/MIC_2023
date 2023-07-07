#!/bin/bash

# Get the bashrc information for conda.
source ~/.bashrc

# Activate the quast conda environment.
conda activate quast_env

# The number of cpu threads for quast.
num_threads=14

# A minimum length of 1500 bps.
minimum_length=1500

# The list of sample ids.
list_file="KGHS_pilot_subset_4_sample_list.txt"

# The input directory to the filtered metagenome files.
input_dir="filtered_metagenomes"

# Internal Field Separator (IFS) used when using cat in a job array so that lines are separated by newlines instead of separated based on spaces.
IFS=$'\n'

for sample_id in $(cat KGHS_pilot_subset_4_sample_list.txt);
do

    # The sample name output directory.
    sample_dir="${input_dir}/${sample_id}"

    # The filtered metagenome assembly fasta file.
    metagenome_fasta_file="${sample_dir}/${sample_id}_min${minimum_length}.fasta"

    # The output directory to write the quast files.
    quast_output_dir="${input_dir}/${sample_id}/quast"

    # Run the quast command to generate the filtered metagenome assembly evaluation metrics file.
    echo "quast.py --output-dir ${quast_output_dir} --threads ${num_threads} ${metagenome_fasta_file}"
    quast.py --output-dir ${quast_output_dir} --threads ${num_threads} ${metagenome_fasta_file}

done

