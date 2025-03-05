#!/bin/bash
#SBATCH -J 1_launch_amplicon_sequence_selected_extractor.py.sh
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL


ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_swarm


python species_comparison_venn.py --taxo_target_gene ../1_EcoPCR/uniq_taxo.txt --taxo_amplicon_gene uniq_taxo.txt --good_discrimination_target_gene ../1_EcoPCR/uniq_taxo_good_discriminated.txt --good_discrimination_amplicon uniq_taxo_good_discriminated.txt --output output_image.png
