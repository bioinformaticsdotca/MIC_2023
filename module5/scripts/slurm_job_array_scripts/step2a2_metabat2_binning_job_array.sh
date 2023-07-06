#!/bin/bash
#SBATCH --job-name="run_metabat2_binning_job_array"
#SBATCH --partition=synergy,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --time=7-00:00:00
#SBATCH --mem=60G
#SBATCH --array=1-4%4
#SBATCH --output=run_metabat2_binning_job_array.%A_%a.out
#SBATCH --error=run_metabat2_binning_job_array.%A_%a.err

# Get the bashrc information for conda.
source ~/.bashrc

# The number of cpu threads for metabat2.
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

# The input directory to write the filtered files.
input_dir="filtered_metagenomes"

# The sample name input directory.
sample_dir="${input_dir}/${sample_id}"

# The filtered metagenome assembly fasta file.
metagenome_fasta_file="${sample_dir}/${sample_id}_min${minimum_length}.fasta"

# The bam file with reads mapped to the filtered metagenome assembly file.
metagenome_bam_file="${sample_dir}/${sample_id}.bam"

# The initial binning directory.
binning_dir="initial_binning"

# The metabat2 depth file.
metabat2_depth_file="${binning_dir}/${sample_id}/${sample_id}_metabat2_depth.txt"

# The metabat2 bin outfile prefix.
metabat2_bin_outfile_prefix="${binning_dir}/${sample_id}/working_dir/metabat2/${sample_id}_bin"

# The metabat2 working directory.
metabat2_working_dir="${binning_dir}/${sample_id}/working_dir/metabat2"
mkdir -p ${metabat2_working_dir}

# The metabat2 bin directory.
metabat2_bin_dir="${binning_dir}/${sample_id}/metabat2"
mkdir -p ${metabat2_bin_dir}

# Activate the metabat2 conda environment.
conda activate metabat2_env

# Run the jgi_summarize_bam_contig_depths command in the metabat2 suite.
echo "jgi_summarize_bam_contig_depths --outputDepth ${metabat2_depth_file} ${metagenome_bam_file}"
jgi_summarize_bam_contig_depths --outputDepth ${metabat2_depth_file} ${metagenome_bam_file}

# Run the metabat2 command to generate initial metabat2 bins.
echo "metabat2 -i ${metagenome_fasta_file} -a ${metabat2_depth_file} -o ${metabat2_bin_outfile_prefix} -m ${minimum_length} -t ${num_threads} --unbinned"
metabat2 -i ${metagenome_fasta_file} -a ${metabat2_depth_file} -o ${metabat2_bin_outfile_prefix} -m ${minimum_length} -t ${num_threads} --unbinned

cp ${metabat2_bin_outfile_prefix}.[0-9]*.fa ${metabat2_bin_dir}

