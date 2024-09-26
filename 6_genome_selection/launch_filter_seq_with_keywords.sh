#!/bin/bash
#SBATCH -J launch_filter_seq_with_keywords
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load devel/python/Python-3.11.1

python filter_seq_with_keywords.py -s /PATH/BD_TaxonMarker/BD_ecoPCR/name_seq_with_taxo.txt -i [included_keywords] -e [excluded_keywords]
