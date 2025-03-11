echo "Starting concatenating and Sorting Results"
combined_results_file="combined_results.tsv"
sorted_results_file="sorted_results.tsv"
filtered_results_file="filtered_data.tsv"
header_file="header.tsv"

cat $(grep -L '^$' result_stat_primers/*_primer_couple.tsv) > "$combined_results_file"


head -n 1 "$combined_results_file" > "$header_file"
sed -i '1!{/^OG_ID/d;}' "$combined_results_file"
# Sorted by column 43 (score), descending order
sort -t$'\t' -k43,43nr "$combined_results_file" > sorted_data.tsv
# We retrieve the 3 best values from column 43
top_scores=$(awk -F'\t' '{print $43}' sorted_data.tsv | sort -nr | uniq | head -n 3)
pattern=$(echo $top_scores | tr ' ' '|')
# Filters out the rows corresponding to these score values
awk -v pattern="$pattern" -F'\t' 'NR==1; NR>1 && $43 ~ pattern' sorted_data.tsv > "$filtered_results_file"

# add header
cat "$header_file" "$filtered_results_file" > "$sorted_results_file"

rm sorted_data.tsv "$combined_results_file" "$header_file" "$filtered_results_file"

echo "completed"



