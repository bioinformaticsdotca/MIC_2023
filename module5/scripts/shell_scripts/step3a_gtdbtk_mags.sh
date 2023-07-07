#!/bin/bash

# Get the bashrc information for conda.
source ~/.bashrc

# The number of cpu threads for metawrap bin refinement module..
num_threads=14

# The completeness threshold
completeness_thresh=50

# The contamination threshold
contamination_thresh=10

# The path to the gtdbtk database.
gtdbtk_data_path="/bulk/IMCshared_bulk/kevin/module5/software_dir/release207_v2"

# The list of sample ids.
list_file="KGHS_pilot_subset_4_sample_list.txt"

# The bin refinement directory for refined bins using metawrap.
bin_refinement_dir="bin_refinement"

# Internal Field Separator (IFS) used when using cat in a job array so that lines are separated by newlines instead of separated based on spaces.
IFS=$'\n'

for sample_id in $(cat $list_file);
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

    # Activate the gtdbtk conda environment.
    conda activate gtdbtk_env

    # Set the path for the gtdbtk database.
    echo "export GTDBTK_DATA_PATH=\"${gtdbtk_data_path}\""
    export GTDBTK_DATA_PATH="${gtdbtk_data_path}"

    # Run gtdbtk on the refined bins from metawrap using a for loop.
    for bin_file in $(ls ${metawrap_bin_refinement_dir} | grep "\.fa")
    do
        echo $bin_file;
        filename=$(basename $bin_file ".fa")
        bin_dir="${sample_refined_bins_dir}/${sample_id}_${filename}"
        mkdir -p $bin_dir

        # The gtdbtk bin output file directory.
        gtdbtk_bin_dir="${bin_dir}/gtdbtk"
        mkdir -p $gtdbtk_bin_dir

        echo "cp ${metawrap_bin_refinement_dir}/${bin_file} ${gtdbtk_bin_dir}/${sample_id}_${bin_file}"
        cp ${metawrap_bin_refinement_dir}/${bin_file} ${gtdbtk_bin_dir}/${sample_id}_${bin_file}
        
        echo "gtdbtk classify_wf --genome_dir ${gtdbtk_bin_dir} --extension \"fa\" --cpus ${num_threads} --out_dir ${gtdbtk_bin_dir}"
        gtdbtk classify_wf --genome_dir ${gtdbtk_bin_dir} --extension "fa" --cpus ${num_threads} --out_dir ${gtdbtk_bin_dir}

        
    done
done
