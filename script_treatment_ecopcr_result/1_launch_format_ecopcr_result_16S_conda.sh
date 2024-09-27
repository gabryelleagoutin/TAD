#!/bin/bash
#SBATCH -J 1_launch_format_ecopcr_result_16S 
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL


ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_swarm

python $PATH_TAXONMARKER/script_treatment_ecopcr_result/format_ecopcr_result.py -o . -t name_seq_16S.txt 16S.ecopcr -vsearch output_vsearch_cluster.txt -taxo name_seq_16S.txt
