#!/bin/bash
#SBATCH --partition=synergy
#SBATCH --nodelist="sm1"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=7-00:00:00
#SBATCH --mem=600G
#SBATCH --output=run_install_DRAM_database.%A.out
#SBATCH --error=run_install_DRAM_database.%A.err

# Get the conda environment path at the start of the shell script.
source ~/.bashrc

# Activate the DRAM conda environment.
conda activate dram_env

# The path to the DRAM database directory you want to use.
dram_database="/bulk/IMCshared_bulk/shared/dbs/DRAM"

# Precomputing the DRAM database. Since we do not have access to the amino acid FASTA file downloaded from KEGG (kegg.pep) because we need a KEGG subscription. We made the database without the use of the KEGG database.
DRAM-setup.py prepare_databases --threads 32 --output_dir ${dram_database}




