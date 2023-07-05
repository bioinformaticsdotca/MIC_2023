#!/bin/bash
#SBATCH --job-name="run_mobileog_db_mags_job_array"
#SBATCH --partition=synergy,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=7-00:00:00
#SBATCH --mem=60G
#SBATCH --array=1-4%4
#SBATCH --output=run_mobileog_db_mags_job_array.%A_%a.out
#SBATCH --error=run_mobileog_db_mags_job_array.%A_%a.err

# Get the bashrc information for conda.
source ~/.bashrc

# The number of cpu threads for metawrap bin refinement module..
num_threads=4

# The completeness threshold
completeness_thresh=50

# The contamination threshold
contamination_thresh=10

mobileog_db_diamond_infile="/bulk/IMCshared_bulk/shared/dbs/mobileOG-db/mobileOG-db-beatrix-1.6.dmnd"

mobileog_db_metadata_infile="/bulk/IMCshared_bulk/shared/dbs/mobileOG-db/mobileOG-db-beatrix-1.6-All.csv"

# Number of Diamond Alignments to Report
kvalue=15

# Maximum E-score
escore=1e-20

# Percent of query coverage.
queryscore=80

# Percent of Identical Matches of samples
pidentvalue=30

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

# Activate the mobileog_db conda environment.	 
conda activate mobileogdb_env

# Run mobileog_db on the refined bins from metawrap using a for loop.
for bin_file in $(ls ${metawrap_bin_refinement_dir} | grep "\.fa")
do 
	echo $bin_file;
	filename=$(basename $bin_file ".fa")
	bin_dir="${sample_refined_bins_dir}/${sample_id}_${filename}"
	mkdir -p $bin_dir

	# The mobileog_db bin output file directory.
	mobileog_db_bin_dir="${bin_dir}/mobileog_db"
	mkdir -p $mobileog_db_bin_dir

	# The prokka protein fasta file.
	prokka_protein_file="${bin_dir}/prokka/${sample_id}_${filename}.faa"

	# Run the diamond command.
	echo "diamond blastp -q ${prokka_protein_file} --db ${mobileog_db_diamond_infile} --outfmt 6 stitle qseqid pident bitscore slen evalue qlen sstart send qstart qend -o "${mobileog_db_bin_dir}/${sample_id}_${filename}.tsv" --threads ${num_threads} -k ${kvalue} -e ${escore} --query-cover ${queryscore} --id ${pidentvalue}"
	diamond blastp -q ${prokka_protein_file} --db ${mobileog_db_diamond_infile} --outfmt 6 stitle qseqid pident bitscore slen evalue qlen sstart send qstart qend -o "${mobileog_db_bin_dir}/${sample_id}_${filename}.tsv" --threads ${num_threads} -k ${kvalue} -e ${escore} --query-cover ${queryscore} --id ${pidentvalue}

	echo "python software_dir/mobileOG-db/mobileOG-pl/mobileOGs-pl-kyanite.py --o ${mobileog_db_bin_dir}/annotation --i \"${mobileog_db_bin_dir}/${sample_id}_${filename}.tsv\" -m ${mobileog_db_metadata_infile}"
	python software_dir/mobileOG-db/mobileOG-pl/mobileOGs-pl-kyanite.py --o ${mobileog_db_bin_dir}/annotation --i "${mobileog_db_bin_dir}/${sample_id}_${filename}.tsv" -m ${mobileog_db_metadata_infile}

done

