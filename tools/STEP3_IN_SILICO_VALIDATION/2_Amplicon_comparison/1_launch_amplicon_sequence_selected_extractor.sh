#!/bin/bash
#SBATCH -J 1_launch_16S_sequence_selected_extractor.py.sh
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL



ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_swarm

python amplicon_sequence_selected_extractor.py --genome_file ../1_EcoPCR/all_genome_name.txt --tsv_file $PATH_BD_[AMPLICON]/[AMPLICON]_with_taxonomy.tsv --fasta_file $PATH_BD_[AMPLICON]/[AMPLICON].fna --output_file filtered_[AMPLICON].fna
