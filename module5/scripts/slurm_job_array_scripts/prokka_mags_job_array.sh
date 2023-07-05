#!/bin/bash
#SBATCH --job-name="run_prokka_mags_job_array"
#SBATCH --partition=synergy,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=7-00:00:00
#SBATCH --mem=60G
#SBATCH --array=1-4%4
#SBATCH --output=run_prokka_mags_job_array.%A_%a.out
#SBATCH --error=run_prokka_mags_job_array.%A_%a.err

# Get the bashrc information for conda.
source ~/.bashrc

# The number of cpu threads for metawrap bin refinement module..
num_threads=4

# The completeness threshold
completeness_thresh=50

# The contamination threshold
contamination_thresh=10

# The list of sample ids.
list_file="KGHS_pilot_subset_4_sample_list.txt"

# Internal Field Separator (IFS) used when using cat in a job array so that lines are separated by newlines instead of separated based on spaces.
IFS=$'\n' 

# Making an bash array based on the list file.
array=($(<$list_file))

# The sample id entry.
sample_id=${array[$SLURM_ARRAY_TASK_ID-1]}

# The bin refinement directory for refined bins using metawrap.
bin_refinement_dir="bin_refinement"

# The sample bin refinement directory for refined bins for each sample.
sample_bin_refinement_dir="${bin_refinement_dir}/${sample_id}"

# The refined bins directory.
refined_bins_dir="refined_bins"
mkdir $refined_bins_dir

# The refined bins directory for each sample.
sample_refined_bins_dir="${refined_bins_dir}/${sample_id}"
mkdir $sample_refined_bins_dir

# The metawrap filtered refined bins directory.
metawrap_bin_refinement_dir="${sample_bin_refinement_dir}/metawrap_${completeness_thresh}_${contamination_thresh}_bins"

# Activate the prokka conda environment.	 
conda activate prokka_env

# Run prokka on the refined bins from metawrap using a for loop.
for bin_file in $(ls ${metawrap_bin_refinement_dir} | grep "\.fa")
do 
	echo $bin_file;
	filename=$(basename $bin_file ".fa")
	bin_dir="${sample_refined_bins_dir}/${sample_id}_${filename}"
	mkdir -p $bin_dir

	# The prokka bin output file directory.
	prokka_bin_dir="${bin_dir}/prokka"
	mkdir -p $prokka_bin_dir

	echo "cp ${metawrap_bin_refinement_dir}/${bin_file} ${prokka_bin_dir}/${sample_id}_${bin_file}"
	cp ${metawrap_bin_refinement_dir}/${bin_file} ${prokka_bin_dir}/${sample_id}_${bin_file}

	refined_bin_file="${metawrap_bin_refinement_dir}/${bin_file}"

	prefix="${sample_id}_${filename}"
	
	echo "prokka --metagenome --outdir ${prokka_bin_dir} --prefix ${prefix} ${refined_bin_file} --cpus ${num_threads} --rfam 1 --force"
	prokka --metagenome --outdir ${prokka_bin_dir} --prefix ${prefix} ${refined_bin_file} --cpus ${num_threads} --rfam 1 --force
	
done

