#!/bin/bash
#SBATCH -J launch_degeprime_trim.
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL
module load devel/python/Python-3.11.1

python process_primers_stat.py -i degeprime_result/concatenated_* -og ../3_fasta_recovery/updated_OG_selected_1578.tab -o result_stat_primers -nm 80 -tm_max 65 -tm_min 54
