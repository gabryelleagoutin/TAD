#!/bin/bash
#SBATCH -J 2_launch_swarm.sh
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load devel/python/Python-3.11.1
module load bioinfo/swarm/3.1.3

python /PATH/TaxonMarker/script_treatment_ecopcr_result/Launch_swarm.py -f all_modified.fna -s fichier_swarm.txt -o cluster.txt  -t 4 -a 1 -d 1


