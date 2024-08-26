#!/bin/bash
#SBATCH -J 3_fasta_recovery
#SBATCH -p workq
#SBATCH --mem=30G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load devel/python/Python-3.11.1 

python fasta_recovery.py ../2_search_taxid_and_monocopy_calculation/OG_[TAXID]_selected.tab
