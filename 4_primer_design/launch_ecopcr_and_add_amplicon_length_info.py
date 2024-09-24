#!/usr/bin/env python

import argparse
import os

__author__ = 'Gabryelle Agoutin - INRAE'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'gabryelle.agoutin@inrae.fr'
__status__ = 'prod'

#script written in python 2.7 to be compatible with EcoPCR in the conda environment


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
    min_sizes = []
    max_sizes = []

    # Read the TSV file line by line
    with open(tsv_file, 'r') as f:
        headers = f.readline().strip().split('\t')
        rows = [line.strip().split('\t') for line in f]

    # Find indices of relevant columns
    og_id_idx = headers.index('OG_ID')
    primer_a_idx = headers.index('Primer_A')
    reverse_complement_b_idx = headers.index('Reverse_Complement_B')

    # Iterate through each row of the TSV file
    for index, row in enumerate(rows):
        og_id = row[og_id_idx]
        primer_a = row[primer_a_idx]
        reverse_complement_b = row[reverse_complement_b_idx]
        
        # Find the appropriate database path
        db_path = find_database_path(og_id, db_paths)
        if db_path:
            # Construct the output file name
            output_file = "{0}_{1}.ecopcr".format(og_id, index + 1)
            command = "ecoPCR -d {0}/{1} {2} {3} > {4}".format(db_path, og_id, primer_a, reverse_complement_b, output_file)
            print("Running: {0}".format(command))
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
            print("Database for {0} not found. Skipping...".format(og_id))
            min_sizes.append(None)
            max_sizes.append(None)

    # Write updated TSV with new columns
    output_tsv_file = tsv_file.replace('.tsv', '_updated.tsv')
    with open(output_tsv_file, 'w') as f_out:
        f_out.write('\t'.join(headers + ['min_size_amplicon', 'max_size_amplicon']) + '\n')
        for i, row in enumerate(rows):
            min_size = str(min_sizes[i]) if min_sizes[i] is not None else ''
            max_size = str(max_sizes[i]) if max_sizes[i] is not None else ''
            f_out.write('\t'.join(row + [min_size, max_size]) + '\n')

    print("Updated TSV file saved as: {0}".format(output_tsv_file))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ecoPCR for each primer pair and analyze the results.",
    epilog='python script.py sorted_results.tsv /path/to/db1 /path/to/db2 or /path/to/db*')
    parser.add_argument("tsv_file", help="Path to the TSV file containing primer information.")
    parser.add_argument("db_paths", nargs='+', help="Paths to the formatted ecoPCR databases.")
    
    args = parser.parse_args()
    
    run_ecopcr(args.tsv_file, args.db_paths)
