#!/bin/bash
#SBATCH --job-name="run_quast_mags_job_array"
#SBATCH --partition=synergy,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --time=7-00:00:00
#SBATCH --mem=60G
#SBATCH --array=1-4%4
#SBATCH --output=run_quast_mags_job_array.%A_%a.out
#SBATCH --error=run_quast_mags_job_array.%A_%a.err

# Get the bashrc information for conda.
source ~/.bashrc

# The number of cpu threads for metawrap bin refinement module..
num_threads=14

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

# Activate the quast conda environment.	 
conda activate quast_env

# Run quast on the refined bins from metawrap using a for loop.
for bin_file in $(ls ${metawrap_bin_refinement_dir} | grep "\.fa")
do 
	echo $bin_file;
	filename=$(basename $bin_file ".fa")
	bin_dir="${sample_refined_bins_dir}/${sample_id}_${filename}"
	mkdir -p $bin_dir

	# The quast bin output file directory.
	quast_bin_dir="${bin_dir}/quast"
	mkdir -p $quast_bin_dir

	echo "cp ${metawrap_bin_refinement_dir}/${bin_file} ${quast_bin_dir}/${sample_id}_${bin_file}"
	cp ${metawrap_bin_refinement_dir}/${bin_file} ${quast_bin_dir}/${sample_id}_${bin_file}

	refined_bin_file="${metawrap_bin_refinement_dir}/${bin_file}"
	
	# Run the quast command to generate the filtered metagenome assembled genome (MAG) evaluation metrics file.
	echo "quast.py --output-dir ${quast_bin_dir} --threads ${num_threads} ${refined_bin_file}"
	quast.py --output-dir ${quast_bin_dir} --threads ${num_threads} ${refined_bin_file}
	
done

