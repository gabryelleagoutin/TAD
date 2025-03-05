#!/bin/bash
#SBATCH -J launch_ecoPCR_16S.sh
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL


ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_EcoPCR

ecoPCR -d ecoPCR_db/[AMPLICON] [PRIMER_FORWARD PRIMER_REVERSE] > [AMPLICON].ecopcr

#ex 16S: ecoPCR -d ecoPCR_db/16S CCTACGGGNGGCWGCAG GACTACHVGGGTATCTAATCC > 16S.ecopcr
