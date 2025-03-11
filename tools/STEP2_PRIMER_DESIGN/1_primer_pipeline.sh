#!/bin/bash

#############################################
# Loading module / activating conda 
#############################################

ml devel/Miniconda/Miniconda3
source $MY_CONDA_PATH
conda activate TaxonMarker_swarm

###############################################################################
# Script: full_pipeline.sh
# Objective: Complete pipeline (alignment, trimming, job generation and launch,
#            DegePrime, stats, linking, and concatenation)
#
# Available parameters:
#   --fasta_path <path>           : Path to the directory containing the .fa files
#                                   (mandatory, no default value).
#   --degeneracies_values <int>   : List of degeneracy values to be used.
#                                   Example: 48 96 
#                                   This is a mandatory parameter.
#
#   --min_primer_length <int>     : Minimum primer length (default: 14).
#   --max_primer_length <int>     : Maximum primer length (default: 24).
#
#   --tab_og_updated <path>       : File for the -og parameter in process_primers_stat.py.
#                                   This is the updated OG file after downloading the FASTA files.
#                                   (default is empty, so it is mandatory if Step 6 is used).
#
#   --number_matching_max <int>   : Parameter -nm for process_primers_stat.py (default: 80).
#                                   It represents the number of sequences that match this primer 
#                                    divided by the total number of sequences * 100. This is the minimum we expect.
#
#   --tm_max <int>                : Parameter -tm_max for process_primers_stat.py (default: 65).
#   --tm_min <int>                : Parameter -tm_min for process_primers_stat.py (default: 54).
#
#   --amplicon_min_size <int>     : Minimum amplicon size (default: 150).
#   --amplicon_max_size <int>     : Maximum amplicon size (default: 590).
#
#   -h / --help                   : Displays this help message.
#
###############################################################################

# default values
FASTA_PATH=""
DEGENERACIES_VALUES_LIST=()
MIN_PRIMER_LENGTH=14
MAX_PRIMER_LENGTH=24
TAB_OG_UPDATED=""
NUMBER_MATCHING_MAX=80
TM_MAX=65
TM_MIN=54
AMPLICON_MIN_SIZE=150
AMPLICON_MAX_SIZE=590


#help

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  --fasta_path <path>            Path to the directory containing the .fa files (mandatory, no default value)."
    echo "  --degeneracies_values <int>    List of degeneracy values (e.g., 48 96). Mandatory."
    echo "  --min_primer_length <int>      Minimum primer length. Default: 14."
    echo "  --max_primer_length <int>      Maximum primer length. Default: 24."
    echo "  --tab_og_updated <path>        File for the -og parameter in process_primers_stat.py. This is the updated OG file after downloading the FASTA files. (default is empty, so it is mandatory if Step 6 is used)."
    echo "  --number_matching_max <int>    Parameter -nm (default: 80).is the number of sequences that match this primer divided by the total number of sequences * 100. This is the minimum we expect."
    echo "  --tm_max <int>                 Parameter -tm_max (default: 65)."
    echo "  --tm_min <int>                 Parameter -tm_min (default: 54)."
    echo "  --amplicon_min_size <int>      Minimum amplicon size (default: 150)."
    echo "  --amplicon_max_size <int>      Maximum amplicon size (default: 590)."
    echo "  -h, --help                     Display this message and exit."
    echo ""
    exit 1
}


# 2) Parsing args
while [[ $# -gt 0 ]]; do
    case $1 in
        --fasta_path)
            FASTA_PATH="$2"
            shift; shift
            ;;
        --degeneracies_values)
            shift
            while [[ $# -gt 0 && ! $1 =~ ^-- ]]; do
                DEGENERACIES_VALUES_LIST+=("$1")
                shift
            done
            ;;
        --min_primer_length)
            MIN_PRIMER_LENGTH="$2"
            shift; shift
            ;;
        --max_primer_length)
            MAX_PRIMER_LENGTH="$2"
            shift; shift
            ;;
        --tab_og_updated)
            TAB_OG_UPDATED="$2"
            shift; shift
            ;;
        --number_matching_max)
            NUMBER_MATCHING_MAX="$2"
            shift; shift
            ;;
        --tm_max)
            TM_MAX="$2"
            shift; shift
            ;;
        --tm_min)
            TM_MIN="$2"
            shift; shift
            ;;
        --amplicon_min_size)
            AMPLICON_MIN_SIZE="$2"
            shift; shift
            ;;
        --amplicon_max_size)
            AMPLICON_MAX_SIZE="$2"
            shift; shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option : $1"
            usage
            ;;
    esac
done

# mandatory args
if [ -z "$FASTA_PATH" ]; then
    echo "Error: you must specify --fasta_path"
    usage
fi

if [ ${#DEGENERACIES_VALUES_LIST[@]} -eq 0 ]; then
    echo "Error: you must specify --degeneracies_values (one or more values)"
    usage
fi

#############################################
# Step 1: Alignment (clustalOmega)
#############################################

echo "Starting Step 1: Alignment"
mkdir -p alignment/

for fichier in "$FASTA_PATH"/*.fa; do
    nom_fichier=$(basename "$fichier")
    echo "Aligning $fichier"
    clustalo -i "$fichier" -o "alignment/$nom_fichier" --threads=4 --force
done
echo "Step 1 completed"

#############################################
# Step 2: Trimming the alignments
#############################################

echo "Starting Step 2: Trimming alignments"
mkdir -p alignment
for fichier in alignment/*.fa; do
    nom_fichier=$(basename "$fichier" .fa)
    echo "Trimming $fichier"
    perl DEGEPRIME/TrimAlignment.pl -i "$fichier" -min 0.9 -o "alignment/trimmed_$nom_fichier.fna"
done
echo "Step 2 completed"

mkdir -p degeprime_result/

#############################################
# Step 3: Creating Sarray files for DegePrime
#############################################

echo "Starting Step 3: Creating Sarray files"

s_array_base='degeprime_multiple_params'

fasta_files='alignment/*.fa'
result_dir='degeprime_result/'
aln_dir='alignment/'

degeneracies_list=("${DEGENERACIES_VALUES_LIST[@]}")
lengths_list=($(seq $MIN_PRIMER_LENGTH 1 $MAX_PRIMER_LENGTH))

degeneracies="${degeneracies_list[*]}"
lengths="${lengths_list[*]}"

line_limit=2450
line_count=0
file_count=1
s_array_file="${s_array_base}_${file_count}.sarray"

echo "# $fasta_files" > "$s_array_file"
echo "# degeneracies $degeneracies" >> "$s_array_file"
echo "# lengths $lengths" >> "$s_array_file"
echo "# result: $result_dir" >> "$s_array_file"

for fasta_file in $fasta_files
do
  cog=$(basename "$fasta_file" .fa)
  trimmed_alignment="$aln_dir/trimmed_${cog}.fna"

  for d in "${degeneracies_list[@]}"
  do
    for l in "${lengths_list[@]}"
    do
      output_fl=$result_dir/${cog}_d${d}_l${l}.tsv

      echo "perl DEGEPRIME/DegePrime.pl -i $trimmed_alignment -d $d -l $l -o $output_fl" >> "$s_array_file"
      ((line_count++))

      # the maximum number of lines for sarray is 2500 lines
      if [ "$line_count" -ge "$line_limit" ]; then
        line_count=0
        ((file_count++))
        s_array_file="${s_array_base}_${file_count}.sarray"
        echo "# $fasta_files" > "$s_array_file"
        echo "# degeneracies $degeneracies" >> "$s_array_file"
        echo "# lengths $lengths" >> "$s_array_file"
        echo "# result: $result_dir" >> "$s_array_file"
      fi
    done
  done
done

echo "Step 3 completed"

#############################################
# Step 4: Launch Sarray for all files generated
#############################################

echo "Starting Step 4: Launching Sarray for all files"

for sarray_file in degeprime_multiple_params_*.sarray; do
    echo "Submitting Sarray file: $sarray_file"
    job_id=$(sarray --mem=200G "$sarray_file" | awk '{print $NF}')
    echo "Submitted batch job $job_id for $sarray_file"

    #  We wait for the job to finish before moving on to the next one
    echo "Waiting for Sarray job $job_id to complete..."
    while squeue -j "$job_id" &> /dev/null; do
        echo "Sarray job for job ID $job_id still running..."
        sleep 60
    done
    echo "Sarray job for job ID $job_id completed."
done

echo "Step 4 completed"

#############################################
# Step 5: Concatenate results
#############################################

echo "Starting Step 5: Concatenating results"
directory="degeprime_result"
ls -l "$directory"

if [ $(ls -1 "$directory"/*.tsv 2>/dev/null | wc -l) -eq 0 ]; then
    echo "No .tsv files found in $directory"
    exit 1
fi

prefixes=$(ls -1 "$directory" | grep -E ".*\.tsv$" | cut -d'_' -f1 | sort | uniq)
echo "Prefixes: $prefixes"

if [ -z "$prefixes" ]; then
    echo "Error: No prefixes found"
    exit 1
fi


for prefix in $prefixes; do
    echo "Processing prefix: $prefix"
    files=$(ls -1 "$directory"/"$prefix"* 2>/dev/null | grep -E ".*\.tsv$")
    echo "Files for prefix $prefix: $files"
    if [ ! -z "$files" ]; then
        echo "Concatenating files for prefix $prefix"
        concatenated_file="$directory/concatenated_$prefix.tsv"
        cat $(echo "$files" | head -n1) > "$concatenated_file"
        for file in $(echo "$files" | tail -n+2); do
            tail -n +2 "$file" >> "$concatenated_file"
        done
    else
        echo "No files found for prefix $prefix"
    fi
done

echo "Step 5 completed"

#############################################
# Step 6: Process Primers Stat
#############################################

echo "Starting Step 6: Processing Primers Stat"
mkdir -p result_stat_primers

if [ -z "$TAB_OG_UPDATED" ]; then
    echo "Warning: --tab_og_updated not specified, the python command below will fail if necessary."
fi

sbatch --wait <<-EOF
#!/bin/bash
#SBATCH -J process_primers_stats
#SBATCH -p unlimitq
#SBATCH --mem=50G
#SBATCH -c 3
#SBATCH --mail-type=BEGIN,END,FAIL

module load devel/python/Python-3.11.1
python process_primers_stat.py \
    -i degeprime_result/concatenated_* \
    -og $TAB_OG_UPDATED \
    -o result_stat_primers \
    -nm $NUMBER_MATCHING_MAX \
    -tm_max $TM_MAX \
    -tm_min $TM_MIN
EOF

echo "Step 6 completed"

#############################################
# Step 7: Creating the Sarray for coupling
#############################################

echo "Starting Step 7: Creating Sarray file for Coupling"
couple_s_array_file='tab_couple.sarray'
rm -f "$couple_s_array_file"

for file in result_stat_primers/concatenated_*; do
    echo "python couple_primer.py -i $file -f alignment/ --amplicon_min_size $AMPLICON_MIN_SIZE --amplicon_max_size $AMPLICON_MAX_SIZE" >> "$couple_s_array_file"
done

echo "Step 7 completed"

#############################################
# Step 8: Launch Sarray for coupling
#############################################

echo "Starting Step 8: Launching Sarray for Coupling"
job_id=$(sarray --mem=200G "$couple_s_array_file" | awk '{print $NF}')
echo "Submitted batch job $job_id"

echo "Waiting for Sarray jobs to complete..."
while squeue -j $job_id &> /dev/null; do
    echo "Sarray jobs still running..."
    sleep 60
done
echo "Sarray jobs completed"

#############################################
# Step 9: Concatenate and sort the final results
#############################################

echo "Starting Step 9: Concatenating and Sorting Results"
combined_results_file="combined_results.tsv"
sorted_results_file="sorted_results.tsv"
filtered_results_file="filtered_data.tsv"
header_file="header.tsv"

cat $(grep -L '^$' result_stat_primers/*_primer_couple.tsv) > "$combined_results_file"

head -n 1 "$combined_results_file" > "$header_file"
sed -i '1!{/^OG_ID/d;}' "$combined_results_file"
sort -t$'\t' -k43,43nr "$combined_results_file" > sorted_data.tsv
top_scores=$(awk -F'\t' '{print $43}' sorted_data.tsv | sort -nr | uniq | head -n 3)
pattern=$(echo $top_scores | tr ' ' '|')
awk -v pattern="$pattern" -F'\t' 'NR==1; NR>1 && $43 ~ pattern' sorted_data.tsv > "$filtered_results_file"
cat "$header_file" "$filtered_results_file" > "$sorted_results_file"

rm sorted_data.tsv "$combined_results_file" "$header_file" "$filtered_results_file"

echo "Step 9 completed"

echo "Pipeline completed successfully!"

