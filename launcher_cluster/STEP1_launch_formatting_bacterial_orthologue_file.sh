#!/bin/bash
#SBATCH -J formatting_bacterial_orthologue_file
#SBATCH -p unlimitq
#SBATCH --mem=100G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load devel/python/Python-3.7.9

python formatting_bacterial_orthologue_file_final.py -o ../Orthodb/odb11v0_OG2genes.tab -b only_line_bacteria.txt -u uniq_og_ids.txt -s ../Orthodb/odb11v0_level2species.tab -g ../Orthodb/odb11v0_OGs.tab -f Bacterial_OG.tab
