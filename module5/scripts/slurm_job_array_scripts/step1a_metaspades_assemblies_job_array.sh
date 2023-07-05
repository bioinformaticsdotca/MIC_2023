#!/bin/bash
#SBATCH --job-name="run_metaspades_assemblies_job_array"
#SBATCH --partition=synergy,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=7-00:00:00
#SBATCH --mem=220G
#SBATCH --array=1-4%4
#SBATCH --output=run_metaspades_assemblies_job_array.%A_%a.out
#SBATCH --error=run_metaspades_assemblies_job_array.%A_%a.err

# Get the bashrc information for conda.
source ~/.bashrc

# Activate the spades conda environment to run metaspades.
conda activate spades_env

# The number of cpu threads for metaspades.
num_threads=28

# The RAM memory in GB for metaspades.
memory_in_gb=220

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

# The path metaspades metagenome assembly sample output directory.
output_dir="metagenome_assemblies/${sample_id}"
mkdir -p $output_dir

# Run the metaspades.py program using the read1 and read2 fastq files 
echo "metaspades.py -k 21,33,55,77,99,127 -t $num_threads -m $memory_in_gb -o $output_dir -1 $fastq_read1 -2 $fastq_read2"
metaspades.py -k 21,33,55,77,99,127 -t $num_threads -m $memory_in_gb -o $output_dir -1 $fastq_read1 -2 $fastq_read2

