#!/usr/bin/env python

import argparse

__author__ = 'Gabryelle Agoutin - INRAE'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'gabryelle.agoutin@inrae.fr'
__status__ = 'prod'

def main(name_seq_file, taxo_file, output_file):
    # Read the taxonomy file and store the information in a dictionary
    taxo_dict = {}
    with open(taxo_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t', 1)
            taxid = parts[0]
            taxonomy = parts[1]
            taxo_dict[taxid] = taxonomy

    # Read the sequence file, add the corresponding taxonomy information
    with open(name_seq_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                parts = line.strip().split('taxid=')
                seq_info = parts[0]
                taxid = parts[1][:-1]  # Remove the semicolon
                taxonomy = taxo_dict.get(taxid, "No taxonomy found")
                outfile.write(f"{seq_info}taxid={taxid}; {taxonomy}\n")
            else:
                outfile.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Complete a sequence file with taxonomy information.')
    parser.add_argument('name_seq_file', type=str, help='File containing sequences without taxonomy')
    parser.add_argument('taxo_file', type=str, help='File containing taxonomy information')
    parser.add_argument('output_file', type=str, help='Output file with completed taxonomy')

    args = parser.parse_args()

    main(args.name_seq_file, args.taxo_file, args.output_file)

