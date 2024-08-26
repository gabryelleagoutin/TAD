#!/bin/bash
s_array_file='ecopcr_commands.sarray'

ecopcr_assembly_dir='/PATH/BD_ecoPCR/ecoPCR_db/assemblies/*.adx'

ecopcr_outdir='primers'
mkdir -p $ecopcr_outdir

COG='[COG_name]'
amorce_f='[forward_primer]'
amorce_r='[reverse_primer]'


for assembly in $ecopcr_assembly_dir
    do
    name=$(basename "$assembly" | sed 's/\.adx$//')

    echo "module load bioinfo/ecoPCR/1.0.1;ecoPCR -d $(echo "$assembly" | sed 's/\.adx$//')  $amorce_f $amorce_r  > primers/${name}.ecopcr;" >> $s_array_file
    done
