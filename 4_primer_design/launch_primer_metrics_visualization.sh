#!/bin/bash
#SBATCH -J launch_primer_metrics_visualization
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_swarm


python primer_metrics_visualization.py -i sorted_results_updated.tsv -o primer_metrics_visualization.html
