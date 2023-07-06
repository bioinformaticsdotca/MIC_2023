#!/bin/bash
#SBATCH --job-name="run_maxbin2_binning_job_array"
#SBATCH --partition=synergy,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --time=7-00:00:00
#SBATCH --mem=60G
#SBATCH --array=1-4%4
#SBATCH --output=run_maxbin2_binning_job_array.%A_%a.out
#SBATCH --error=run_maxbin2_binning_job_array.%A_%a.err

# Get the bashrc information for conda.
source ~/.bashrc

# The number of cpu threads for maxbin2.
num_threads=14

# A minimum length of 1500 bps.
minimum_length=1500

# The maxbin2 marker set number.
maxbin2_markers=107

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

# The maxbin2 depth file.
maxbin2_depth_file="${binning_dir}/${sample_id}/${sample_id}_maxbin2_depth.txt"

# The maxbin2 abundance file.
maxbin2_abund_file="${binning_dir}/${sample_id}/${sample_id}_maxbin2_abund.txt"

# The maxbin2 abundance list file.
maxbin2_abund_list_file="${binning_dir}/${sample_id}/${sample_id}_maxbin2_abund_list.txt"

# The maxbin2 bin output file prefix to name the bins.
maxbin2_bin_outfile_prefix="${binning_dir}/${sample_id}/working_dir/maxbin2/${sample_id}_bin"

# The maxbin2 working directory.
maxbin2_working_dir="${binning_dir}/${sample_id}/working_dir/maxbin2"

# Create the maxbin2 working directory
mkdir -p ${maxbin2_working_dir}

# The maxbin2 bin directory.
maxbin2_bin_dir="${binning_dir}/${sample_id}/maxbin2"

# Create the maxbin2 bin directory
mkdir -p ${maxbin2_bin_dir}

# Activate the metabat2 conda environment.        
conda activate metabat2_env

# Run the jgi_summarize_bam_contig_depths command in the metabat2 suite.
echo "jgi_summarize_bam_contig_depths --outputDepth ${maxbin2_depth_file} --noIntraDepthVariance ${metagenome_bam_file}"
jgi_summarize_bam_contig_depths --outputDepth ${maxbin2_depth_file} --noIntraDepthVariance ${metagenome_bam_file}

# Use tail to skip the header line using tail -n+2 and grab the first and third column and write to the maxbin2 abundance file. 
echo "tail -n+2 ${maxbin2_depth_file} | cut -f1,3 > ${maxbin2_abund_file}"
tail -n+2 ${maxbin2_depth_file} | cut -f1,3 > ${maxbin2_abund_file}

# Add the maxbin2 abundance file path to the maxbin2 abundance list file for maxbin2. 
echo "echo "${maxbin2_abund_file}" > ${maxbin2_abund_list_file}"
echo "${maxbin2_abund_file}" > ${maxbin2_abund_list_file}

# Activate the maxbin2 conda environment.	 
conda activate maxbin2_env

# Run the maxbin2 command.
echo "perl software_dir/MaxBin-2.2.7/run_MaxBin.pl -contig ${metagenome_fasta_file} -markerset ${maxbin2_markers} -thread ${num_threads} -min_contig_length ${minimum_length} -out ${maxbin2_bin_outfile_prefix} -abund_list ${maxbin2_abund_list_file}"
perl software_dir/MaxBin-2.2.7/run_MaxBin.pl -contig ${metagenome_fasta_file} -markerset ${maxbin2_markers} -thread ${num_threads} -min_contig_length ${minimum_length} -out ${maxbin2_bin_outfile_prefix} -abund_list ${maxbin2_abund_list_file}

# Rename and copy the maxbin2 bins to the maxbin2 bin directory.
for bin_file in $(ls ${maxbin2_working_dir} | grep "\.fasta");
do 
	echo $bin_file;
 
    	# Get the filename using basename to remove the ".fasta" extension.
	filename=$(basename $bin_file '.fasta');
    
    	# Get the base number for the bin.
    	bin_num=$(echo $filename | sed -r "s/${sample_id}_bin\.0+//g");
    	echo $bin_num;
	
    	# Make the new filename using the sample id and the bin number.
    	new_filename="${sample_id}_bin.${bin_num}.fa"
    	echo $new_filename;
    
    	# Copy the file to the maxbin2 folder and rename the file.
	cp ${maxbin2_working_dir}/$bin_file ${maxbin2_bin_dir}/$new_filename;
    
done

         


