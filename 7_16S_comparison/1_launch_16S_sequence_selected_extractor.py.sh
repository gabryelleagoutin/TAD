#!/bin/bash
#SBATCH -J launch_format_ecoPCR_result
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load devel/python/Python-3.11.1

python 16S_sequence_selected_extractor.py --genome_file ../6_genome_selection/genome_name_selected.txt --tsv_file /home/gagoutin/work/BD_TaxonMarker/BD_16S/16S_with_taxonomy.tsv --fasta_file /home/gagoutin/work/BD_TaxonMarker/BD_16S/16S.fna --output_file filtered_16S.fna
