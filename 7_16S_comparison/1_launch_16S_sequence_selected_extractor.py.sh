#!/bin/bash
#SBATCH -J 1_launch_16S_sequence_selected_extractor.py.sh
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load devel/python/Python-3.11.1

python 16S_sequence_selected_extractor.py --genome_file ../6_genome_selection/genome_name_selected.txt --tsv_file /PATH/BD_TaxonMarker/BD_16S/16S_with_taxonomy.tsv --fasta_file /PATH/BD_TaxonMarker/BD_16S/16S.fna --output_file filtered_16S.fna
