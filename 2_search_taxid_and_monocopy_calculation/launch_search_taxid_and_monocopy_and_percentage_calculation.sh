#!/bin/bash
#SBATCH -J search_taxid_and_monocopy_calculation
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 20
#SBATCH --mail-type=BEGIN,END,FAIL

module load devel/python/Python-3.7.9

python search_taxid_and_monocopy_and_percentage_calculation.py -i ../1_formatting_file_bacterian_OG/save_results/Bacterial_OG.tab -f ../Orthodb/odb11v0_species.tab -l ../Orthodb/odb11v0_level2species.tab -s 1578 -o OG_1760.tab
