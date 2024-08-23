#!/bin/bash
s_array_file='ecopcr_commands.sarray'

ecopcr_assembly_dir='/home/gagoutin/work/BD_ecoPCR/ecoPCR_db/assemblies/*.adx'

ecopcr_outdir='primers'
COG='58151at1578'
amorce_f='AAYATGGGKGTBGGNGAYG'
amorce_r='TTCWGGGAARAYBARYTGT'


for assembly in $ecopcr_assembly_dir
    do
    name=$(basename "$assembly" | sed 's/\.adx$//')

    echo "module load bioinfo/ecoPCR/1.0.1;ecoPCR -d $(echo "$assembly" | sed 's/\.adx$//')  $amorce_f $amorce_r  > primers/${name}.ecopcr;" >> $s_array_file
    done
