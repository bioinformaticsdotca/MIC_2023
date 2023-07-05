#!/bin/bash
#SBATCH --job-name="run_quast_metagenome_metrics_job_array"
#SBATCH --partition=synergy,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --time=7-00:00:00
#SBATCH --mem=60G
#SBATCH --array=1-4%4
#SBATCH --output=run_quast_metagenome_metrics_job_array.%A_%a.out
#SBATCH --error=run_quast_metagenome_metrics_job_array.%A_%a.err

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

# Internal Field Separator (IFS) used when using cat in a job array so that lines are separated by newlines instead of separated based on spaces.
IFS=$'\n' 

# Making an bash array based on the list file.
array=($(<$list_file))

# The sample id entry.
sample_id=${array[$SLURM_ARRAY_TASK_ID-1]}

# A minimum length of 1500 bps. 
minimum_length=1500

# The input directory to the filtered metagenome files.
input_dir="filtered_metagenomes"

# The sample name output directory.
sample_dir="${input_dir}/${sample_id}"

# The filtered metagenome assembly fasta file.
metagenome_fasta_file="${sample_dir}/${sample_id}_min${minimum_length}.fasta"

# The output directory to write the quast files.
quast_output_dir="${input_dir}/${sample_id}/quast"

# Run the quast command to generate the filtered metagenome assembly evaluation metrics file.
echo "quast.py --output-dir ${quast_output_dir} --threads ${num_threads} ${metagenome_fasta_file}"
quast.py --output-dir ${quast_output_dir} --threads ${num_threads} ${metagenome_fasta_file}


