#!/bin/bash
#SBATCH -J 3_launch_stat_swarm_16S
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load devel/python/Python-3.11.1

python /PATH/TaxonMarker/script_treatment_ecopcr_result/stats_report_taxo2.py -c cluster.txt -o stats.txt -l rejected_clusters.txt -r cluster_corrected.txt

