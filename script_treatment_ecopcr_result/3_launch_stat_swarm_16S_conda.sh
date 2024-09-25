#!/bin/bash
#SBATCH -J 3_launch_stat_swarm_16S
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL


ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_swarm

python $PATH_TAXONMARKER/script_treatment_ecopcr_result/stats_report_taxo2.py -c cluster.txt -o stats.txt -l rejected_clusters.txt -r cluster_corrected.txt

