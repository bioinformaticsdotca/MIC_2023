#!/bin/bash
#SBATCH --job-name="run_map_reads_filtered_assembly_job_array"
#SBATCH --partition=synergy,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=7-00:00:00
#SBATCH --mem=60G
#SBATCH --array=1-4%4
#SBATCH --output=run_map_reads_filtered_assembly_job_array.%A_%a.out
#SBATCH --error=run_map_reads_filtered_assembly_job_array.%A_%a.err

# Get the bashrc information for conda.
source ~/.bashrc

# Activate the bwa conda environment.
conda activate bwa_env

# The number of cpu threads for BWA and samtools..
num_threads=4

# A minimum length of 1500 bps.
minimum_length=1500

# The list of sample ids.
list_file="KGHS_pilot_subset_4_sample_list.txt"

# Internal Field Separator (IFS) used when using cat in a job array so that lines are separated by newlines instead of separated based on spaces.
IFS=$'\n' 

# Making an bash array based on the list file.
array=($(<$list_file))

# The sample id entry.
sample_id=${array[$SLURM_ARRAY_TASK_ID-1]}

# The path fastq input file directory.
input_dir="cleaned_fastq_files"

# The read 1 fastq file for each sample.
fastq_read1="${input_dir}/${sample_id}_1.fastq"

# The read 2 fastq file for each sample.
fastq_read2="${input_dir}/${sample_id}_2.fastq"

# A minimum length of 1500 bps. 
minimum_length=1500

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

