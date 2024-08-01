#!/bin/bash
#SBATCH -J launch_alignment
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load bioinfo/ClustalOmega/1.2.4

mkdir -p alignment/

for fichier in ../3_fasta_recovery/*.fa; do
    nom_fichier=$(basename "$fichier")
    clustalo -i "$fichier" -o "alignment/$nom_fichier" --threads=4
done


