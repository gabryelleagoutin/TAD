#!/bin/bash
for file in result_stat_primers/concatenated_*; do echo "module load devel/python/Python-3.11.1;python couple_primer.py -i $file -f alignment/ --amplicon_min_size 150 --amplicon_max_size 590" >> tab_couple.sarray;done
