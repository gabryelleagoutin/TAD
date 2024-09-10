#!/usr/bin/env python

import argparse
import os
import pandas as pd

__author__ = 'Gabryelle Agoutin - INRAE'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'gabryelle.agoutin@inrae.fr'
__status__ = 'prod'



def find_database_path(og_id, db_paths):
    """
    Find the database path that matches the given OG_ID.

    Parameters:
    og_id (str): The OG_ID to search for.
    db_paths (list of str): List of paths to the formatted ecoPCR databases.

    Returns:
    str: Path to the database that matches the OG_ID, or None if not found.
    """
    for path in db_paths:
        if og_id in os.path.basename(path):
            return path
    return None

def extract_amplicon_sizes(ecopcr_file):
    """
    Extract amplicon sizes from the ecoPCR output file.

    Parameters:
    ecopcr_file (str): Path to the ecoPCR output file.

    Returns:
    list of int: List of amplicon sizes extracted from the file.
    """
    sizes = []
    with open(ecopcr_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.split(' | ')
            if len(parts) > 19:
                try:
                    size = int(parts[19].strip())
                    sizes.append(size)
                except ValueError:
                    continue
    return sizes

def run_ecopcr(tsv_file, db_paths):
    """
    Run ecoPCR for each line in the TSV file and analyze the output files.

    Parameters:
    tsv_file (str): Path to the TSV file containing primer information.
    db_paths (list of str): List of paths to the formatted ecoPCR databases.
    """
    # Load the TSV file into a DataFrame
    df = pd.read_csv(tsv_file, sep='\t')
    
    min_sizes = []
    max_sizes = []

    # Iterate through each row of the DataFrame
    for index, row in df.iterrows():
        og_id = row['OG_ID']
        primer_a = row['Primer_A']
        reverse_complement_b = row['Reverse_Complement_B']
        
        # Find the appropriate database path
        db_path = find_database_path(og_id, db_paths)
        if db_path:
            # Construct the output file name
            output_file = f"{og_id}_{index + 1}.ecopcr"
            command = f"ecoPCR -d {db_path}/{og_id} {primer_a} {reverse_complement_b} > {output_file}"
            print(f"Running: {command}")
            os.system(command)

            # Analyze the output file
            sizes = extract_amplicon_sizes(output_file)
            if sizes:
                min_sizes.append(int(min(sizes)))
                max_sizes.append(int(max(sizes)))
            else:
                min_sizes.append(None)
                max_sizes.append(None)
        else:
            print(f"Database for {og_id} not found. Skipping...")
            min_sizes.append(None)
            max_sizes.append(None)

    # Add the min and max sizes to the DataFrame
    df['min_size_amplicon'] = min_sizes
    df['max_size_amplicon'] = max_sizes

    # Convert columns to integer type (to remove '.0')
    df['min_size_amplicon'] = df['min_size_amplicon'].fillna(0).astype(int)
    df['max_size_amplicon'] = df['max_size_amplicon'].fillna(0).astype(int)

    # Save the updated DataFrame to a new TSV file
    output_tsv_file = tsv_file.replace('.tsv', '_updated.tsv')
    df.to_csv(output_tsv_file, sep='\t', index=False)
    print(f"Updated TSV file saved as: {output_tsv_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ecoPCR for each primer pair and analyze the results.",
    epilog='python script.py sorted_results.tsv /path/to/db1 /path/to/db2 or /path/to/db*')
    parser.add_argument("tsv_file", help="Path to the TSV file containing primer information.")
    parser.add_argument("db_paths", nargs='+', help="Paths to the formatted ecoPCR databases.")
    
    args = parser.parse_args()
    
    run_ecopcr(args.tsv_file, args.db_paths)
