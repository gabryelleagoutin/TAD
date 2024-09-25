#!/bin/bash
#SBATCH -J 3_launch_stat_swarm
#SBATCH -p workq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL


ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_swarm

python $PATH_TAXONMARKER/script_treatment_ecopcr_result/stats_report_taxo.py -c analysis_1_GGNGAAMGDACNCGTGAAG-TCNGAHARTTCRTCCATMCC/cluster.txt -o analysis_1_GGNGAAMGDACNCGTGAAG-TCNGAHARTTCRTCCATMCC/stats.txt -i f__Lactobacillaceae -e g__Leuconostoc g__Pediococcus g__Unassigned sp. g__Weissella g__Dellaglioa g__Periweissella g__Fructobacillus -l analysis_1_GGNGAAMGDACNCGTGAAG-TCNGAHARTTCRTCCATMCC/rejected_clusters.txt -r analysis_1_GGNGAAMGDACNCGTGAAG-TCNGAHARTTCRTCCATMCC/cluster_corrected.txt

