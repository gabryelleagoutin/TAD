#!/bin/bash

# Define the path of the directory containing the files
directory="degeprime_result"

# Recover file prefixes
prefixes=$(ls -1 "$directory" | grep -E ".*\.tsv$" | cut -d'_' -f1 | sort | uniq)

# Loop on each prefix

for prefix in $prefixes; do
    # Concatenate all files with the same prefix
    files=$(ls -1 "$directory"/"$prefix"* 2>/dev/null | grep -E ".*\.tsv$")
    if [ ! -z "$files" ]; then
        concatenated_file="$directory/concatenated_$prefix.tsv"
        cat $(echo "$files" | head -n1) > "$concatenated_file"
        for file in $(echo "$files" | tail -n+2); do
            tail -n +2 "$file" >> "$concatenated_file"
        done
    fi
done
