#!/bin/bash

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

    # Activate the dram conda environment.	 
    conda activate dram_env

    # Run dram on the refined bins from metawrap using a for loop.
    for bin_file in $(ls ${metawrap_bin_refinement_dir} | grep "\.fa")
    do 
        echo $bin_file;
        filename=$(basename $bin_file ".fa")
        bin_dir="${sample_refined_bins_dir}/${sample_id}_${filename}"
        mkdir -p $bin_dir

        # The dram bin output file directory.
        dram_bin_dir="${bin_dir}/dram"
        #mkdir -p $dram_bin_dir

        refined_bin_file="${metawrap_bin_refinement_dir}/${bin_file}"

        echo "DRAM.py annotate -i ${refined_bin_file} -o ${dram_bin_dir}"
        DRAM.py annotate -i ${refined_bin_file} -o ${dram_bin_dir}

        echo "DRAM.py distill -i ${dram_bin_dir}/annotations.tsv -o ${dram_bin_dir}/genome_summaries --trna_path ${dram_bin_dir}/trnas.tsv --rrna_path ${dram_bin_dir}/rrnas.tsv"
        DRAM.py distill -i ${dram_bin_dir}/annotations.tsv -o ${dram_bin_dir}/genome_summaries --trna_path ${dram_bin_dir}/trnas.tsv --rrna_path ${dram_bin_dir}/rrnas.tsv
        
    done

done
