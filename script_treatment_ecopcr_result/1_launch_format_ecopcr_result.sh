#!/bin/bash
#SBATCH -J 1_launch_format_ecopcr_result
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load devel/python/Python-3.11.1
ml bioinfo/VSEARCH/2.22.1

python /PATH/TaxonMarker/script_treatment_ecopcr_result/format_ecopcr_result.py -o . -t /PATH/BD_TaxonMarker/BD_ecoPCR/name_seq_with_taxo.txt primers/*ecopcr
