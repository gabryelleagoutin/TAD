#!/usr/bin/env python3

"""
Merge assemblies and format them to be processed by obiconvert to build an ecoPCR database.

see also: format_sequences.py
"""

# Metadata
__author__ = 'Mainguy Jean - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2020 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.1'
__email__ = 'support.bioinfo.genotoul@inra.fr'
__status__ = 'dev'

import gzip
import os
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys
import csv
from Bio import SeqIO
import re
from collections import defaultdict

def tsv_to_dict_of_dicts(file, key_field):
    """
    Take a tsv with header and parse into dict of dict
    Uses
    * the specified field as key
    * the line turned into a dict as value .
    """
    dict_of_list_of_dict = defaultdict(list)
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for l in reader:
            dict_of_list_of_dict[l[key_field]].append(l)
        return dict(dict_of_list_of_dict)

def get_genomic_fna_fl(assembly_root_dir, assembly_accession):
    """
    Retrieve genomic fna file for a given assembly accession.
    """

    preassembly_path = assembly_accession.replace('_', '').split('.')[0]
    assembly_path = '/'.join([preassembly_path[index: index + 3]
                              for index in range(0, len(preassembly_path), 3)])

    assembly_dir = os.path.join(assembly_root_dir, assembly_path)

    if not os.path.isdir(assembly_dir):
        raise ValueError(f'{assembly_dir} does not exists...')

    genomic_fna_fls = []
    for f in os.listdir(assembly_dir):
        if f.startswith(assembly_accession) and f.endswith("_genomic.fna.gz"):
            genomic_fna_fls.append(f)

    if len(genomic_fna_fls) != 1:
        raise ValueError(f'There is more than one file than ends with _genomic.fna.gz in {assembly_dir} dir')
    else:
        genomic_fna_fl = genomic_fna_fls.pop()
    return os.path.join(assembly_dir, genomic_fna_fl)

def parse_arguments():
    parser = ArgumentParser(description="Merge assemblies and format them to be processed by obiconvert to build an ecoPCR database.",
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument("assembly_list_to_keep", type=str, default=None,
                        help="")
    parser.add_argument("--assembly_selection", type=str,
                        help="", default=None)
    parser.add_argument("--assembly_root_dir", type=str,
                        help="", default='ncbi_data/refseq_assemblies/')

    parser.add_argument("-o", "--output", type=str, default=None,
                        help="")
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")

    args = parser.parse_args()
    return args

def main():

    args = parse_arguments()

    if args.verbose:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)
        logging.info('Mode verbose ON')

    else:
        logging.basicConfig(format="%(levelname)s: %(message)s")

    assembly_root_dir = args.assembly_root_dir
    assembly_list_to_keep = args.assembly_list_to_keep
    assemblies_with_taxonomy = args.assembly_selection

    if args.output:
        merged_assemblies_file = args.output
    else:
        merged_assemblies_file = f"{os.path.basename(assembly_list_to_keep).split('.')[0]}_merged.fna"

    assemblies_id_mapper_fl = '.'.join(merged_assemblies_file.split('.')[:-1]) + '_map.ids'

    assembly_to_keep = set()

    with open(assembly_list_to_keep) as fl:
        assembly_to_keep = {a.rstrip() for a in fl}
    logging.info(f'{len(assembly_to_keep)} assemblies to format and merged')

    with open(assemblies_with_taxonomy) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        assembly_to_taxid = {d['assembly_accession']: d['taxid']
                             for d in reader if d['assembly_accession'] in assembly_to_keep}

    logging.info('Retrieving fna files of assemblies')
    assembly_to_fna_fl = {a: get_genomic_fna_fl(
        assembly_root_dir, a) for a in assembly_to_keep}

    logging.info(f'Formatting fasta header, merging sequences and writing them in {merged_assemblies_file}')
    with open(merged_assemblies_file, 'w') as out, open(assemblies_id_mapper_fl, 'w') as map_id_wrt:

        for assembly_name, fna_fl in assembly_to_fna_fl.items():
            proper_open = gzip.open if fna_fl.endswith('.gz') else open
            with proper_open(fna_fl, "rt") as fl:

                for i, seq in enumerate(SeqIO.parse(fl, "fasta")):
                    taxid = assembly_to_taxid[assembly_name]
                    original_id = seq.id
                    new_id = f'{assembly_name}|{i+1:04d}'
                    if len(str(i+1)) > 4:
                        logging.critical(f'there are at least {i+1} contigs in assembly {assembly_name}. 4 trailing 0 was not enough')
                    seq.description = f'taxid={taxid}; {seq.description}'

                    seq.id = new_id

                    SeqIO.write(seq, out, "fasta")
                    map_id_wrt.write(f'{new_id}\t{original_id}\n')

if __name__ == '__main__':
    main()
