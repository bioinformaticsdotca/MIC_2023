#!/bin/bash

# Get the conda environment path at the start of the shell script.
source ~/.bashrc

# Activate the DRAM conda environment.
conda activate dram_env

# The path to the DRAM database directory you want to use.
dram_database="/bulk/IMCshared_bulk/shared/dbs/DRAM"

# Precomputing the DRAM database. Since we do not have access to the amino acid FASTA file downloaded from KEGG (kegg.pep) because we need a KEGG subscription. We made the database without the use of the KEGG database.
DRAM-setup.py prepare_databases --threads 32 --output_dir ${dram_database}




