#!/bin/bash

# Get the bashrc information for conda.
source ~/.bashrc

# Activate the bwa conda environment.
conda activate bwa_env

# The number of cpu threads for BWA and samtools.
num_threads=14

# A minimum length of 1500 bps.
minimum_length=1500

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

    # The output directory to write the filtered files.
    output_dir="filtered_metagenomes"

    # The sample name output directory.
    sample_dir="${output_dir}/${sample_id}"

    # The filtered metagenome assembly fasta file.
    metagenome_fasta_file="${sample_dir}/${sample_id}_min${minimum_length}.fasta"

    # The sam file with reads mapped to the filtered metagenome assembly file.
    metagenome_sam_file="${sample_dir}/${sample_id}.sam"

    # Index the metagenome fasta file as the reference to map reads.
    echo "bwa index ${metagenome_fasta_file}"
    bwa index ${metagenome_fasta_file}

    # Map the reads to the filtered metagenome using bwa mem for binning.
    echo "bwa mem -t ${num_threads} ${metagenome_fasta_file} ${fastq_read1} ${fastq_read2} > ${metagenome_sam_file}"
    bwa mem -t ${num_threads} ${metagenome_fasta_file} ${fastq_read1} ${fastq_read2} > ${metagenome_sam_file}

    # The filtered metagenome bam file.
    metagenome_bam_file="${sample_dir}/${sample_id}.bam"

    # Activate the samtools conda environment.
    conda activate samtools_env

    # Run samtools to sort the sam file and convert to bam.
    echo "samtools sort -@ ${num_threads} -O BAM -o ${metagenome_bam_file} ${metagenome_sam_file}"
    samtools sort -@ ${num_threads} -O BAM -o ${metagenome_bam_file} ${metagenome_sam_file}

done
