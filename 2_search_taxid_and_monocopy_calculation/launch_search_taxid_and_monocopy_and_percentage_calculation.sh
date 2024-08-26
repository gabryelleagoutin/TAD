#!/bin/bash
#SBATCH -J 2_search_taxid_and_monocopy_calculation
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 20
#SBATCH --mail-type=BEGIN,END,FAIL

module load devel/python/Python-3.11.1

python search_taxid_and_monocopy_and_percentage_calculation.py -i ../1_formatting_file_bacterian_OG/save_results/Bacterial_OG.tab -f PATH/Orthodb/odb11v0_species.tab -l PATH/Orthodb/odb11v0_level2species.tab -s [TAXID] -o OG_[TAXID].tab
