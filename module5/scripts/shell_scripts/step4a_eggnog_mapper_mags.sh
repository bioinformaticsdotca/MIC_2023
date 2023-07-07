#!/bin/bash

# Get the bashrc information for conda.
source ~/.bashrc

# The number of cpu threads for metawrap bin refinement module..
num_threads=14

# The completeness threshold
completeness_thresh=50

# The contamination threshold
contamination_thresh=10

# The eggnog_mapper_db folder path.
eggnog_mapper_db="/bulk/IMCshared_bulk/shared/shared_software/eggnog-mapper/data"

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

    # Activate the eggnog_mapper conda environment.
    conda activate eggnog_mapper_env

    # Run eggnog_mapper on the refined bins from metawrap using a for loop.
    for bin_file in $(ls ${metawrap_bin_refinement_dir} | grep "\.fa")
    do
        echo $bin_file;
        filename=$(basename $bin_file ".fa")
        bin_dir="${sample_refined_bins_dir}/${sample_id}_${filename}"
        mkdir -p $bin_dir

        # The eggnog_mapper bin output file directory.
        eggnog_mapper_bin_dir="${bin_dir}/eggnog_mapper"
        mkdir -p $eggnog_mapper_bin_dir

        # The prokka protein fasta file.
        prokka_protein_file="${bin_dir}/prokka/${sample_id}_${filename}.faa"

        # The eggnog-mapper output file prefix.
        eggnog_mapper_output_file_prefix="${eggnog_mapper_bin_dir}/${sample_id}_${filename}"

        # Run the run_eggnog_mapper command.
        echo "python /bulk/IMCshared_bulk/shared/shared_software/eggnog-mapper/emapper.py -i ${prokka_protein_file} --itype proteins --cpu ${num_threads} --data_dir ${eggnog_mapper_db} --output ${eggnog_mapper_output_file_prefix}"
        python /bulk/IMCshared_bulk/shared/shared_software/eggnog-mapper/emapper.py -i ${prokka_protein_file} --itype proteins --cpu ${num_threads} --data_dir ${eggnog_mapper_db} --output ${eggnog_mapper_output_file_prefix}
    done
done
