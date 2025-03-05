#!/bin/bash
#SBATCH -J launch_degeprime_trim.
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

for files in alignment/*.fa; do
    file_name=$(basename "$fichier" .fa)
    perl DEGEPRIME/TrimAlignment.pl -i "$files" -min 0.9 -o "alignment/trimmed_$file_name.fna"
done
