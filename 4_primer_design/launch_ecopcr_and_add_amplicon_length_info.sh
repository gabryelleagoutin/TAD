#!/bin/bash
#SBATCH -J launch_ecopcr_and_add_amplicon_length_info
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load devel/python/Python-3.11.1
module load bioinfo/ecoPCR/1.0.1


python launch_ecopcr_and_add_amplicon_length_info.py sorted_results.tsv ecoPCR_db_*
