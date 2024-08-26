#!/bin/bash
#SBATCH -J launch_alignment
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load bioinfo/ClustalOmega/1.2.4

mkdir -p alignment/

for files in ../3_fasta_recovery/*.fa; do
    file_name=$(basename "$fichier")
    clustalo -i "$files" -o "alignment/$file_name" --threads=4
done


