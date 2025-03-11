#!/bin/bash
#SBATCH -J 1_launch_format_ecopcr_result_conda
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_swarm

#GENERAL
python ../script_treatment_ecopcr_result/format_ecopcr_result.py -o analysis_[PRIMER_FOLDER] -t $PATH_BD_ECOPCR/name_seq_with_taxo.txt result/[PRIMER_FOLDER]/*ecopcr

## for example, if there are several files to execute, you can put several commands.
## python ../script_treatment_ecopcr_result/format_ecopcr_result.py -o analysis_1_GGNGAAMGDACNCGTGAAG-TCNGAHARTTCRTCCATMCC -t $PATH_BD_ECOPCR/name_seq_with_taxo.txt result/1_GGNGAAMGDACNCGTGAAG-TCNGAHARTTCRTCCATMCC/*ecopcr
## python ../script_treatment_ecopcr_result/format_ecopcr_result.py -o analysis_2_GNGAAMGDACNCGTGAAGG-TCNGAHARTTCRTCCATMCC -t $PATH_BD_ECOPCR/name_seq_with_taxo.txt result/2_GNGAAMGDACNCGTGAAGG-TCNGAHARTTCRTCCATMCC/*ecopcr
