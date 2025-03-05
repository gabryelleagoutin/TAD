#!/bin/bash
#SBATCH -J 3_fasta_recovery
#SBATCH -p workq
#SBATCH --mem=30G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL


ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_swarm

python fastas_recovery.py ../1_search_taxid_and_monocopy_calculation/OG_[TAXID]_selected.tab
