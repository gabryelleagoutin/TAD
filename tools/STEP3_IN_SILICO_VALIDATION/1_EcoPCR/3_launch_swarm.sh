#!/bin/bash
#SBATCH -J 2_launch_swarm.sh
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL


ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_swarm

python ../script_treatment_ecopcr_result/Launch_swarm.py -f [PRIMER_FOLDER]/all_modified.fna -s [PRIMER_FOLDER]/fichier_swarm.txt -o [PRIMER_FOLDER]/cluster.txt  -t 4 -a 1 -d 1 -vsearch output_vsearch_cluster.txt -taxo $PATH_BD_ECOPCR/name_seq_with_taxo.txt


