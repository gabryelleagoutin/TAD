#!/bin/bash

s_array_file='degeprime_multiple_params.sarray'

fasta_files='alignment/*.fa' 
result_dir='degeprime_result/'
aln_dir='alignment/'

degeneracies='8 12 24 48 96'
lengths='14 15 16 17 18 19 20 21 22 23 24'

echo '#' $fasta_files > $s_array_file
echo '# degeneracies' $degeneracies >> $s_array_file
echo '# lengths' $lengths >> $s_array_file
echo '# result:' $result_dir >> $s_array_file

for fasta_file in $fasta_files
do
  cog=$(basename "$fasta_file" .fa)
  trimmed_alignment="$aln_dir/trimmed_${cog}.fna"
  
  for d in $degeneracies
  do
    for l in $lengths
    do
      output_fl=$result_dir/${cog}_d${d}_l${l}.tsv

      echo "perl DEGEPRIME/DegePrime.pl -i $trimmed_alignment -d $d -l $l -o $output_fl" >> $s_array_file
    done
  done
done

echo sarray -J degeprime -o slurm_array_out/%j_%x.out $s_array_file
