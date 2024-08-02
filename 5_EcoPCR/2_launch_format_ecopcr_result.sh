#!/bin/bash
#SBATCH -J launch_format_ecoPCR_result
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load devel/python/Python-3.11.1
ml bioinfo/VSEARCH/2.22.1

python format_ecopcr_result.py -o . -t name_seq.txt primers/*ecopcr
