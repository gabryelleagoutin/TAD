#!/bin/bash
#SBATCH -J launch_iPCR2FrogsDB
#SBATCH -p unlimitq
#SBATCH --mem=50G

ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_swarm


python iPCRresults_to_valide_file_for_FrogsDB.py --fasta all.fna --taxonomy $PATH_BD_ECOPCR/name_seq_with_taxo.txt  --amplicon_name [AMPLICON_NAME] --creation_date YYYYMMDD
