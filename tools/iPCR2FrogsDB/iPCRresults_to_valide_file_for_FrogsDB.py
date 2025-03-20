#!/usr/bin/env python3

__author__ = 'Gabryelle Agoutin - INRAE & Jean Mainguy - Genoscope'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'gabryelle.agoutin@inrae.fr'
__status__ = 'prod'


import argparse

def parse_taxonomy_file(taxonomy_file):
    """
    Builds a dictionary { ID: taxonomy } from a file with lines like :
        >GCA_000686505.1|0001 taxid=1339247; k__Bacteria;p__...
    With 3 "columns" separated by spaces.

    The 1st column => '>GCA_000686505.1|0001'
    The 2nd column => 'taxid=1339247;'
    The 3rd column => 'k__Bacteria;p__Bacillota;...'

    """
    taxonomy_dict = {}
    with open(taxonomy_file, 'r') as f:
        for line in f:
            line = line.strip()

            if not line or not line.startswith('>'):
                continue
                
            parts = line.split()
            first_part = parts[0][1:]
            seq_id = first_part.split('|')[0]
            if len(parts) > 2:
                taxonomy = ' '.join(parts[2:])
            else:
                taxonomy = ''

            taxonomy_dict[seq_id] = taxonomy

    return taxonomy_dict

def annotate_and_extract(fasta_file, taxonomy_dict, cleaned_fasta, tax_file):
    """
    Reads a fasta file and produces two outputs :
      1) a new FASTA with a cleaned seqID (without "|" parts)
      2) a file listing seqID + taxonmy.
    """
    with open(fasta_file, 'r') as fin, \
         open(cleaned_fasta, 'w') as fout_fasta, \
         open(tax_file, 'w') as fout_tax:

        current_header = None
        current_seq_lines = []

        for line in fin:
            line = line.rstrip('\n')
            if line.startswith(">"):
                if current_header is not None and current_seq_lines:
                    seq_str = "\n".join(current_seq_lines)
                    fout_fasta.write(f"{current_header}\n{seq_str}\n")

                # New sequence
                raw_id = line[1:].split("|")[0]
                current_header = f">{raw_id}"

                # Taxo retrived from the dictionary created earlier
                tax = taxonomy_dict.get(raw_id, "")

                fout_tax.write(f"{raw_id}\t{tax}\n")

                current_seq_lines = []
            else:
                current_seq_lines.append(line)

        # last seq
        if current_header is not None and current_seq_lines:
            seq_str = "\n".join(current_seq_lines)
            fout_fasta.write(f"{current_header}\n{seq_str}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Creates a new FASTA with cleaned headers and a taxonomy file (ID, taxonomy).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--fasta", required=True, help="FULL ecoPCR output fasta file after the ‘format_ecopcr_result’ step without dereplication (ex: all.fna).")
    parser.add_argument("--taxonomy", required=True, help="Taxonomy file (ex: name_seq_with_taxo.txt).")
    parser.add_argument("--amplicon_name", required=True, help="Amplicon name (ex: 16S).")
    parser.add_argument("--creation_date", required=True, help="Creation date YYYYMMDD (ex: 20250320).")

    args = parser.parse_args()


    taxonomy_dict = parse_taxonomy_file(args.taxonomy)

    cleaned_fasta_output = f"{args.amplicon_name}_db_{args.creation_date}.fasta"
    tax_file_output = f"{args.amplicon_name}_db_{args.creation_date}.tax"

    annotate_and_extract(
        fasta_file=args.fasta,
        taxonomy_dict=taxonomy_dict,
        cleaned_fasta=cleaned_fasta_output,
        tax_file=tax_file_output
    )

    print(f"Files generated:\n"
          f"  - FASTA : {cleaned_fasta_output}\n"
          f"  - Taxonomy : {tax_file_output}")


if __name__ == "__main__":
    main()
