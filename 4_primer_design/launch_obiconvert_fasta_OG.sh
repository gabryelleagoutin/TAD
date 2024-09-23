#!/bin/bash

# Load necessary modules
ml devel/Miniconda/Miniconda3
ml bioinfo/OBITools/1.2.11

# Check if at least two arguments are provided (path to ncbi_tax_dump and at least one FASTA file)
if [ "$#" -lt 2 ]; then
  echo "Usage: $0 /path/to/ncbi_tax_dump fasta_file1 [fasta_file2 ...] or all fasta *.fa.
  the best use is to select the OG fastas in the results table. 
  Formatting is a fairly long process, so select the ones you have rather than formatting them all."
  exit 1
fi

# The first argument is the path to the ncbi_tax_dump database
ncbi_tax_dump_path=$1

# All subsequent arguments are FASTA files
shift
fasta_files=("$@")

# Function to run obiconvert with the given FASTA file and ncbi_tax_dump path
run_obiconvert() {
  local fasta_file=$1
  local ncbi_tax_dump_path=$2
  
  # Extract the OG name by removing "_fasta.fa" from the filename
  local og_name=$(basename "$fasta_file" | sed 's/_fasta\.fa//')
  
  # Build the output path for --ecopcrdb-output
  local ecopcrdb_output="ecoPCR_db_${og_name}/${og_name}"
  
  # Create the output directory if it does not exist
  local output_dir=$(dirname "$ecopcrdb_output")
  if [ ! -d "$output_dir" ]; then
    echo "Creating output directory: $output_dir"
    mkdir -p "$output_dir"
  fi
  
  # Build the obiconvert command
  local command="obiconvert --fasta $fasta_file --ecopcrdb-output=$ecopcrdb_output -t $ncbi_tax_dump_path"
  
  # Execute the obiconvert command
  echo "Executing: $command"
  eval $command
  
  # Check if the command executed successfully
  if [ $? -ne 0 ]; then
    echo "Error executing $command"
    exit 1  # Exit script if command fails
  fi
}

# Iterate over each FASTA file and run the obiconvert function
for fasta_file in "${fasta_files[@]}"; do
  run_obiconvert "$fasta_file" "$ncbi_tax_dump_path"
done
