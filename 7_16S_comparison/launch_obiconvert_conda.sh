#!/bin/bash
#SBATCH -J launch_obiconvert
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL


ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_EcoPCR

mkdir ecoPCR_db/
obiconvert --fasta filtered_16S.fna --ecopcrdb-output=ecoPCR_db/16S -t $PATH_NCBI_TAX_DUMB
