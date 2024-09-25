#!/bin/bash
#SBATCH -J 1_launch_16S_sequence_selected_extractor.py.sh
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL



ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_swarm

python 16S_sequence_selected_extractor.py --genome_file $PATH_TAXONMARKER/6_genome_selection/genome_name_selected.txt --tsv_file $PATH_BD_16S/16S_with_taxonomy.tsv --fasta_file $PATH_BD_16S/16S.fna --output_file filtered_16S.fna
