#!/bin/bash

#SBATCH -J full_pipeline
#SBATCH -p unlimitq
#SBATCH --mem=100G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

# Load necessary modules
module load bioinfo/ClustalOmega/1.2.4

# Step 1: Alignment with clustalOmega
echo "Starting Step 1: Alignment"
mkdir -p alignment/

for fichier in ../3_fasta_recovery/*.fa; do
    nom_fichier=$(basename "$fichier")
    echo "Aligning $fichier"
    clustalo -i "$fichier" -o "alignment/$nom_fichier" --threads=4 --force
done
echo "Step 1 completed"

# Step 2: Trimming alignments
echo "Starting Step 2: Trimming alignments"
for fichier in alignment/*.fa; do
    nom_fichier=$(basename "$fichier" .fa)
    echo "Trimming $fichier"
    perl DEGEPRIME/TrimAlignment.pl -i "$fichier" -min 0.9 -o "alignment/trimmed_$nom_fichier.fna"
done
echo "Step 2 completed"

mkdir -p degeprime_result/

# Step 3: Create Sarray file for DegePrime
echo "Starting Step 3: Creating Sarray file"
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

echo "Step 3 completed"

# Step 4: Launch Sarray
echo "Starting Step 4: Launching Sarray"
job_id=$(sarray --mem=200G degeprime_multiple_params.sarray | awk '{print $NF}')
echo "Submitted batch job $job_id"

# Wait for Sarray jobs to finish
echo "Waiting for Sarray jobs to complete..."
while squeue -j $job_id > /dev/null 2>&1; do
    echo "Sarray jobs still running..."
    sleep 60
done
echo "Sarray jobs completed"

# Step 5: Concatenate results
echo "Starting Step 5: Concatenating results"
directory="degeprime_result"
ls -l "$directory"

# Verify files existence before proceeding
if [ $(ls -1 "$directory"/*.tsv 2>/dev/null | wc -l) -eq 0 ]; then
    echo "No .tsv files found in $directory"
    exit 1
fi

# Retrieve file prefixes
prefixes=$(ls -1 "$directory" | grep -E ".*\.tsv$" | cut -d'_' -f1 | sort | uniq)
echo "Prefixes: $prefixes"

# Check if prefixes is empty
if [ -z "$prefixes" ]; then
    echo "Error: No prefixes found"
    exit 1
fi

# Loop through each prefix
for prefix in $prefixes; do
    echo "Processing prefix: $prefix"
    # Concatenate all files with the same prefix
    files=$(ls -1 "$directory"/"$prefix"* 2>/dev/null | grep -E ".*\.tsv$")
    echo "Files for prefix $prefix: $files"
    if [ ! -z "$files" ]; then
        echo "Concatenating files for prefix $prefix"
        concatenated_file="$directory/concatenated_$prefix.tsv"
        # Concatenate all files, removing the header except for the first file
        cat $(echo "$files" | head -n1) > "$concatenated_file"
        for file in $(echo "$files" | tail -n+2); do
            tail -n +2 "$file" >> "$concatenated_file"
        done
    else
        echo "No files found for prefix $prefix"
    fi
done

echo "Step 5 completed"

# Step 6: Process Primers Stat

# -og and parameter need to be changed 
echo "Starting Step 6: Processing Primers Stat"
mkdir -p result_stat_primers
sbatch --wait <<-EOF
#!/bin/bash
#SBATCH -J process_primers_stats
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL
module load devel/python/Python-3.11.1
python process_primers_stat.py -i degeprime_result/concatenated_* -og ../3_fasta_recovery/updated_OG_1760_selected.tab -o result_stat_primers -nm 80 -tm_max 65 -tm_min 54
EOF
echo "Step 6 completed"

# Step 7: Create Sarray file for Coupling
echo "Starting Step 7: Creating Sarray file for Coupling"
couple_s_array_file='tab_couple.sarray'
rm -f $couple_s_array_file
for file in result_stat_primers/concatenated_*; do
    echo "module load devel/python/Python-3.11.1;python couple_primer.py -i $file -f alignment/ --amplicon_min_size 150 --amplicon_max_size 590" >> $couple_s_array_file
done

echo "Step 7 completed"

# Step 8: Launch Sarray for Coupling
echo "Starting Step 8: Launching Sarray for Coupling"
job_id=$(sarray --mem=200G $couple_s_array_file | awk '{print $NF}')
echo "Submitted batch job $job_id"

# Wait for Sarray jobs to finish
echo "Waiting for Sarray jobs to complete..."
while squeue -j $job_id > /dev/null 2>&1; do
    echo "Sarray jobs still running..."
    sleep 60
done
echo "Sarray jobs completed"

# Step 9: Concatenate and Sort Results
echo "Starting Step 9: Concatenating and Sorting Results"
combined_results_file="combined_results.tsv"
sorted_results_file="sorted_results.tsv"
filtered_results_file="filtered_data.tsv"
header_file="header.tsv"

cat $(grep -L '^$' result_stat_primers/*_primer_couple.tsv) > $combined_results_file
head -n 1 $combined_results_file > $header_file
sed -i '1!{/^OG_ID/d;}' $combined_results_file
sort -t$'\t' -k43,43nr $combined_results_file > sorted_data.tsv
top_scores=$(awk -F'\t' '{print $43}' sorted_data.tsv | sort -nr | uniq | head -n 3)
pattern=$(echo $top_scores | tr ' ' '|')
awk -v pattern="$pattern" -F'\t' 'NR==1; NR>1 && $43 ~ pattern' sorted_data.tsv > $filtered_results_file
cat $header_file $filtered_results_file > $sorted_results_file
rm sorted_data.tsv $combined_results_file $header_file

echo "Step 9 completed"

echo "Pipeline completed successfully!"
