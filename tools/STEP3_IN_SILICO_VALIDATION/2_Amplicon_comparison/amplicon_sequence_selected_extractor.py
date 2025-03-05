#!/usr/bin/env python

import argparse
import pandas as pd
import logging

__author__ = 'Gabryelle Agoutin - INRAE'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'gabryelle.agoutin@inrae.fr'
__status__ = 'prod'

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################


def setup_logging(log_file):
    """
    Sets up logging configuration to output messages to the specified log file.

    Parameters:
    log_file (str): Path to the log file where log messages will be saved.
    """
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def read_genome_names(file_path):
    """
    Reads a list of genome names from the specified file.

    Parameters:
    file_path (str): Path to the file containing genome names.

    Returns:
    set: A set of genome names.
    """
    with open(file_path, 'r') as f:
        return set(line.strip() for line in f)

def read_tsv_file(tsv_file, genome_names):
    """
    Reads the TSV file and filters rows based on the provided genome names.

    Parameters:
    tsv_file (str): Path to the TSV file containing taxonomic information.
    genome_names (set): A set of genome names to filter the TSV file.

    Returns:
    dict: A dictionary mapping seqid to a tuple of (taxid, taxonomy).
    """
    df = pd.read_csv(tsv_file, sep='\t')
    filtered_df = df[df['genome_name'].isin(genome_names)]
    return {
        row['seqid']: (row['species_taxid'], row['taxonomy']) for _, row in filtered_df.iterrows()
    }

def filter_fasta(fasta_file, seqid_info_map, output_file, name_file):
    """
    Filters sequences from the FASTA file based on seqid information, writes the result to an output file,
    and writes the full sequence header (including taxonomy) to a separate file.

    Parameters:
    fasta_file (str): Path to the input FASTA file.
    seqid_info_map (dict): A dictionary mapping seqid to a tuple of (taxid, taxonomy).
    output_file (str): Path to the output FASTA file where filtered sequences will be saved.
    name_file (str): Path to the output file where full sequence headers will be saved.
    """
    with open(fasta_file, 'r') as f_in, open(output_file, 'w') as f_out, open(name_file, 'w') as name_out:
        write_seq = False
        for line in f_in:
            if line.startswith('>'):
                parts = line.split('|')
                seqid = parts[0][1:]
                if seqid in seqid_info_map:
                    taxid, taxonomy = seqid_info_map[seqid]
                    # Header for FASTA file (without taxonomy)
                    header_fasta = f'>{seqid}| taxid={taxid};'
                    # Header for name file (with taxonomy)
                    header_name = f'>{seqid}| taxid={taxid}; {taxonomy}'
                    f_out.write(f'{header_fasta}\n')
                    name_out.write(f'{header_name}\n')
                    write_seq = True
                else:
                    write_seq = False
            elif write_seq:
                f_out.write(line)

    logging.info(f"Selection and extraction completed. Filtered sequences have been saved to {output_file}.")
    logging.info(f"Sequence names (including taxonomy) have been saved to {name_file}.")

def extract_unique_taxonomies(name_file, uniq_taxonomy_file):
    """
    Extracts unique taxonomies from the name file and writes them to a unique taxonomy file.

    Parameters:
    name_file (str): Path to the file containing sequence headers with taxonomies.
    uniq_taxonomy_file (str): Path to the output file where unique taxonomies will be saved.
    """
    taxonomies = set()

    with open(name_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Extract the taxonomy part of the header
                parts = line.split('; ')
                if len(parts) > 1:
                    taxonomy = parts[1].strip()
                    taxonomies.add(taxonomy)

    with open(uniq_taxonomy_file, 'w') as f:
        for taxonomy in sorted(taxonomies):
            f.write(f'{taxonomy}\n')

    logging.info(f"Unique taxonomies have been saved to {uniq_taxonomy_file}.")
    logging.info(f"Number of unique taxonomies: {len(taxonomies)}.")

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################


def main():
    """
    Main function to coordinate reading input files, filtering sequences, writing the output files,
    generating a file with full sequence headers, and extracting unique taxonomies.
    """
    parser = argparse.ArgumentParser(
        description='Filter amplicons sequences based on selected genomes and extract unique taxonomies.',
        epilog='Example usage: python script.py --genome_file genome_name_selected.txt --tsv_file 16S_with_taxonomy.tsv --fasta_file 16S.fna --output_file filtered_16S.fna'
    )
    parser.add_argument('--genome_file', required=True, help='Path to the file containing selected genome names.')
    parser.add_argument('--tsv_file', required=True, help='Path to the TSV file containing taxonomic information.')
    parser.add_argument('--fasta_file', required=True, help='Path to the input FASTA file containing amplicon sequences.')
    parser.add_argument('--output_file', required=True, help='Path to the output file for filtered sequences.')
    parser.add_argument('--name_file', default='name_seq_amplicon.txt', help='Path to the output file where sequence names with taxonomy will be saved.')
    parser.add_argument('--log_file', default='process.log', help='Path to the log file where messages will be saved.')

    args = parser.parse_args()

    setup_logging(args.log_file)

    genome_names = read_genome_names(args.genome_file)
    seqid_info_map = read_tsv_file(args.tsv_file, genome_names)

    filter_fasta(args.fasta_file, seqid_info_map, args.output_file, args.name_file)

    uniq_taxonomy_file = 'uniq_taxonomy_amplicon.txt'
    extract_unique_taxonomies(args.name_file, uniq_taxonomy_file)

if __name__ == "__main__":
    main()
