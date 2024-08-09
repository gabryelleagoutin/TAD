#!/bin/bash
#SBATCH -J launch_format_ecoPCR_result
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load bioinfo/ecoPCR/1.0.1

ecoPCR -d ecoPCR_db/16S CCTACGGGNGGCWGCAG GACTACHVGGGTATCTAATCC > 16S.ecopcr
