#!/bin/bash
#SBATCH -J launch_ecopcr_and_add_amplicon_length_info
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL


ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_EcoPCR

python ecopcr_and_add_amplicon_length_info.py ../sorted_results.tsv ecoPCR_db_*
