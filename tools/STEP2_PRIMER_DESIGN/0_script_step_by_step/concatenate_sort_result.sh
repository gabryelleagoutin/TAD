cat $(grep -L '^$' result_stat_primers/*_primer_couple.tsv) > combined_results.tsv
head -n 1 combined_results.tsv > header.tsv
sed -i '1!{/^OG_ID/d;}' combined_results.tsv
sort -t$'\t' -k44,44nr combined_results.tsv > sorted_data.tsv
top_scores=$(awk -F'\t' '{print $44}' sorted_data.tsv | sort -nr | uniq | head -n 3)
pattern=$(echo $top_scores | tr ' ' '|')
awk -v pattern="$pattern" -F'\t' 'NR==1; NR>1 && $44 ~ pattern' sorted_data.tsv > filtered_data.tsv
cat header.tsv filtered_data.tsv > sorted_results.tsv
rm sorted_data.tsv combined_results.tsv header.tsv
