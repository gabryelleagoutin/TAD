#!/bin/bash
#SBATCH -J launch_degeprime_trim.
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

for fichier in alignment/*.fa; do
    nom_fichier=$(basename "$fichier" .fa)
    perl DEGEPRIME/TrimAlignment.pl -i "$fichier" -min 0.9 -o "alignment/trimmed_$nom_fichier.fna"
done
