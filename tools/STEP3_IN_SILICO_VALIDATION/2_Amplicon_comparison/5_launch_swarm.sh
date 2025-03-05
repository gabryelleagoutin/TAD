#!/bin/bash
#SBATCH -J 3_launch_stat_swarm
#SBATCH -p workq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL


ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_swarm

python ../script_treatment_ecopcr_result/stats_report_taxo.py -c [PRIMER_FOLDER]/cluster.txt -o [PRIMER_FOLDER]/stats.txt -i [included_keywords] -e [excluded_keywords] -l [PRIMER_FOLDER]/rejected_clusters.txt -r [PRIMER_FOLDER]/cluster_corrected.txt [--clean_words] [--all_species] [--selected_species]
