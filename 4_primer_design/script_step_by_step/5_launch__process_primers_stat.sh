#!/bin/bash
#SBATCH -J process_primers_stat
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load devel/python/Python-3.11.1

python process_primers_stat.py -i degeprime_result/concatenated_* -og ../3_fasta_recovery/updated_OG_selected_[TAXID].tab -o result_stat_primers -nm [VALUES]  -tm_max [VALUES] -tm_min [VALUES]
