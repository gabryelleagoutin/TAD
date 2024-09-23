#!/bin/bash
#SBATCH -J 1_formatting_file_bacterian_OG
#SBATCH -p unlimitq
#SBATCH --mem=100G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load devel/python/Python-3.11.1

python formatting_bacterial_orthologue_file.py -o PATH/Orthodb/odb11v0_OG2genes.tab -b only_line_bacteria.txt -u uniq_og_ids.txt -s PATH/Orthodb/odb11v0_level2species.tab -g PATH/Orthodb/odb11v0_OGs.tab -f Bacterial_OG.tab
