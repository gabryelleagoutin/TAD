#!/usr/bin/env python

__author__ = 'Gabryelle Agoutin - INRAE & Jean Mainguy - Genoscope'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'gabryelle.agoutin@inrae.fr'
__status__ = 'prod'


import os
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import defaultdict
from Bio.Seq import Seq
from subprocess import run


def get_seq_from_ecopcr_file(ecopcrfile, add_primer_sequences):
    """Get amplicon sequences from an ecopcr file."""
    assembly_counter = defaultdict(int)
    with open(ecopcrfile, encoding="latin-1") as fl:
        for i, l in enumerate(fl):
            if l.startswith('#'):
                continue
            splitted_l = [e.strip() for e in l.split(' | ')]
            if len(splitted_l) < 20:
                logging.critical(f'ecopcr line {i} ({l}) in file {ecopcrfile} is incorrect')
                raise IndexError(f'ecopcr line {i} ({l}) in file {ecopcrfile} is incorrect')
            assembly = splitted_l[0].split('|')[0]
            assembly_counter[assembly] += 1
            f_primer = splitted_l[13]
            r_primer = splitted_l[16]
            if add_primer_sequences:
                r_primer_seq = Seq(r_primer)
                sequence = f_primer + splitted_l[20] + str(r_primer_seq.reverse_complement())
            else:
                sequence = splitted_l[20]
            header = f'{assembly}|seq{assembly_counter[assembly]}'
            yield header, sequence


def write_ecopcr_file_seq_to_fasta(ecopcrfile, output_file, add_primer_sequences):
    """Write fasta sequences from an ecopcr file."""
    logging.info(f'Writing sequences in {output_file}')
    with open(output_file, 'w') as fl:
        for header, seq in get_seq_from_ecopcr_file(ecopcrfile, add_primer_sequences):
            fl.write(f'>{header}\n')
            fl.write(seq + '\n')


def parse_taxonomy_file(taxonomy_file):
    taxonomy_dict = {}
    with open(taxonomy_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                parts = line.split(" ", 1)
                seq_id = parts[0][1:].split("|")[0]
                taxonomy_info = parts[1].strip()
                taxonomy_dict[seq_id] = taxonomy_info
    return taxonomy_dict


def replace_fasta_headers_and_check_sequences(fasta_file, taxonomy_dict, output_file, log_file):
    valid_nucleotides = set("ACTG")
    id_count = defaultdict(int)
    unique_id_set = set()
    seq_id = None
    seq = []
    
    with open(fasta_file, 'r') as fin, open(output_file, 'w') as fout, open(log_file, 'w') as flog:
        for line in fin:
            if line.startswith(">"):
                if seq_id and seq:
                    seq_str = "".join(seq)
                    invalid_chars = set(seq_str) - valid_nucleotides
                    if invalid_chars:
                        invalid_char_str = ", ".join(invalid_chars)
                        flog.write(f"Sequence {seq_id} contains invalid characters: {invalid_char_str}\n")
                        flog.write(f"Sequence: {seq_str}\n\n")
                    else:
                        # Ensure the ID is unique by adding a suffix if necessary
                        unique_id = seq_id
                        while unique_id in unique_id_set:
                            id_count[seq_id] += 1
                            unique_id = f"{seq_id}.{id_count[seq_id]}"
                        unique_id_set.add(unique_id)
                        fout.write(f">{unique_id}| {taxonomy_dict.get(seq_id, '')}\n")
                        fout.write(f"{seq_str}\n")
                seq_id = line[1:].split("|")[0]
                seq = []
            else:
                seq.append(line.strip())
        
        if seq_id and seq:
            seq_str = "".join(seq)
            invalid_chars = set(seq_str) - valid_nucleotides
            if invalid_chars:
                invalid_char_str = ", ".join(invalid_chars)
                flog.write(f"Sequence {seq_id} contains invalid characters: {invalid_char_str}\n")
                flog.write(f"Sequence: {seq_str}\n\n")
            else:
                unique_id = seq_id
                while unique_id in unique_id_set:
                    id_count[seq_id] += 1
                    unique_id = f"{seq_id}.{id_count[seq_id]}"
                unique_id_set.add(unique_id)
                fout.write(f">{unique_id}| {taxonomy_dict.get(seq_id, '')}\n")
                fout.write(f"{seq_str}\n")


def main():
    parser = ArgumentParser(description="Process EcoPCR output and annotate with taxonomy",
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument("ecopcr_files", type=str, nargs='+', help="EcoPCR output files")
    parser.add_argument('--add_primer_sequences', help="Add primer sequences to amplicon sequences", action="store_true")
    parser.add_argument("-o", '--output_dir', help="Output directory", type=str, default='./')
    parser.add_argument("-t", '--taxonomy_file', help="Taxonomy file", type=str, required=True)
    parser.add_argument("-v", "--verbose", help="Increase output verbosity", action="store_true")
    
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)
        logging.info('Verbose mode ON')
    else:
        logging.basicConfig(format="%(levelname)s: %(message)s")

    ecopcrfiles = [f for f in args.ecopcr_files if f.endswith('ecopcr')]
    logging.info(f'{len(ecopcrfiles)} ecopcr files to convert')

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    add_primer_sequences = args.add_primer_sequences

    for ecopcrfile in ecopcrfiles:
        logging.info(f'Processing {ecopcrfile}')
        output_file = os.path.join(output_dir, ''.join(os.path.basename(ecopcrfile).split('.')[:-1]) + '.fna')
        write_ecopcr_file_seq_to_fasta(ecopcrfile, output_file, add_primer_sequences=add_primer_sequences)

    combined_fasta = os.path.join(output_dir, 'all.fna')
    with open(combined_fasta, 'w') as fout:
        for ecopcrfile in ecopcrfiles:
            fasta_file = os.path.join(output_dir, ''.join(os.path.basename(ecopcrfile).split('.')[:-1]) + '.fna')
            with open(fasta_file, 'r') as fin:
                fout.write(fin.read())

    derep_fasta = os.path.join(output_dir, 'derep.fasta')
    run(['vsearch', '--derep_fulllength', combined_fasta, '--sizein', '--sizeout', '--output', derep_fasta])
    taxonomy_dict = parse_taxonomy_file(args.taxonomy_file)
    final_output = os.path.join(output_dir, 'all_modified.fna')
    log_file = os.path.join(output_dir, 'invalid_sequences.log')
    replace_fasta_headers_and_check_sequences(derep_fasta, taxonomy_dict, final_output, log_file)

if __name__ == '__main__':
    main()
