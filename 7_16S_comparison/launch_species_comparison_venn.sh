#!/bin/bash
#SBATCH -J 1_launch_16S_sequence_selected_extractor.py.sh
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL


ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_swarm


python script_name.py --taxo_target_gene path/to/uniq_taxo.txt --taxo_16S_gene path/to/uniq_taxo_16S.txt --good_discrimination_target_gene path/to/uniq_taxo_good_discriminated.txt --good_discrimination_16S path/to/uniq_taxo_good_discriminated_16S.txt --output path/to/output_image.png
