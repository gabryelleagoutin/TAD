#!/bin/bash
#SBATCH -J launch_obiconvert
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

mkdir ecoPCR_db
ml devel/Miniconda/Miniconda3
ml bioinfo/OBITools/1.2.11

mkdir ecoPCR_db/
obiconvert --fasta filtered_16S.fna --ecopcrdb-output=ecoPCR_db/16S -t /PATH/BD_TaxonMarker/ncbi_tax_dumb/2024-15-04
