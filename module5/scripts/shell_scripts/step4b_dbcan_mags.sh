#!/bin/bash

# Get the bashrc information for conda.
source ~/.bashrc

# The number of cpu threads for metawrap bin refinement module..
num_threads=14

# The completeness threshold
completeness_thresh=50

# The contamination threshold
contamination_thresh=10

# The run_dbcan_database folder path.
run_dbcan_database="/bulk/IMCshared_bulk/kevin/module5/software_dir/run_dbcan_database"

# The list of sample ids.
list_file="KGHS_pilot_subset_4_sample_list.txt"

# The bin refinement directory for refined bins using metawrap.
bin_refinement_dir="bin_refinement"

# Internal Field Separator (IFS) used when using cat in a job array so that lines are separated by newlines instead of separated based on spaces.
IFS=$'\n'

for sample_id in $(cat KGHS_pilot_subset_4_sample_list.txt);
do

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

    # Activate the dbcan conda environment.
    conda activate dbcan_env

    # Run dbcan on the refined bins from metawrap using a for loop.
    for bin_file in $(ls ${metawrap_bin_refinement_dir} | grep "\.fa")
    do
        echo $bin_file;
        filename=$(basename $bin_file ".fa")
        bin_dir="${sample_refined_bins_dir}/${sample_id}_${filename}"
        mkdir -p $bin_dir

        # The dbcan bin output file directory.
        dbcan_bin_dir="${bin_dir}/dbcan"
        mkdir -p $dbcan_bin_dir

        # The prokka protein fasta file.
        prokka_protein_file="${bin_dir}/prokka/${sample_id}_${filename}.faa"

        # The refined bin fasta file.
        refined_bin_file="${metawrap_bin_refinement_dir}/${bin_file}"

        # Run the run_dbcan command.
        echo "run_dbcan ${prokka_protein_file} protein --db_dir ${run_dbcan_database} --dia_cpu ${num_threads} --hmm_cpu ${num_threads} --dbcan_thread ${num_threads} --out_dir ${dbcan_bin_dir}"
        run_dbcan ${prokka_protein_file} protein --db_dir ${run_dbcan_database} --dia_cpu ${num_threads} --hmm_cpu ${num_threads} --dbcan_thread ${num_threads} --out_dir ${dbcan_bin_dir}
        
    done
done
