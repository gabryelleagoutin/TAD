#!/bin/bash
#SBATCH -J launch_format_ecoPCR_result
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

mkdir ecoPCR_db
ml devel/Miniconda/Miniconda3
ml bioinfo/OBITools/1.2.11

obiconvert --fasta filtered_16S.fna --ecopcrdb-output=ecoPCR_db/16S -t /home/gagoutin/work/BD_TaxonMarker/ncbi_tax_dumb/2024-15-04
