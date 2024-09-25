#!/bin/bash
#SBATCH -J launch_ecoPCR_16S.sh
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load bioinfo/ecoPCR/1.0.1

ecoPCR -d ecoPCR_db/16S [PRIMER_FOWARD PRIMER_REVERSE] > 16S.ecopcr
