#!/bin/bash
#SBATCH -J 2_search_taxid_and_monocopy_calculation
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 20
#SBATCH --mail-type=BEGIN,END,FAIL


ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_swarm



python search_taxid_and_monocopy_and_percentage_calculation.py -i .../../0_File_bacterian_OG/Bacterial_OG.tab -f $PATH_ORTHODB/Orthodb/odb11v0_species.tab -l $PATH_ORTHODB/Orthodb/odb11v0_level2species.tab -s [TAXID] -o OG_[TAXID].tab
