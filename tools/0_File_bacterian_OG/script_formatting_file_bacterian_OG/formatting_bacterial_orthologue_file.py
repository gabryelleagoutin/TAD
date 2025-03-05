#!/usr/bin/env/ python

__author__ = 'Gabryelle Agoutin - INRAE'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'gabryelle.agoutin@inrae.fr'
__status__ = 'prod'


import sys
import argparse
import subprocess


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

def extract_line_bacteria(OrthoDB_file, bacteria_line_file, species_file):
    '''
    Extracts lines from the OrthoDB_file that match the identifiers
    specified in the bacteria_IDs_file and writes them to the
    bacteria_line_file.

    Parameters:
        OrthoDB_file (str): Path to the OrthoDB file from which lines
            will be extracted.
        bacteria_line_file (str): Path to the output file where matching
            lines will be written.
        species_file (str): Path to the file containing bacterial species identifiers.

    This function reads the identifiers from species_file and
    scans each line in OrthoDB_file. If an identifier in OrthoDB_file
    matches any of the identifiers in species_file, the line is
    written to the bacteria_line_file.
    '''
    command = f"awk '$1==2' {species_file} | cut -f2"
    result = subprocess.check_output(command, shell=True, text=True)

    identifiants = set(result.strip().split('\n'))

    with open(bacteria_line_file, 'w') as blf:
        with open(OrthoDB_file, 'r') as of:
            for line in of:
                parts = line.split()
                if len(parts) >= 2:
                    identifiant = parts[1].split(':')[0]
                    if identifiant in identifiants:
                        blf.write(line)

    print("Finished. The corresponding lines have been written in", bacteria_line_file)


def extract_unique_og_ids(bacteria_line_file, uniq_OG):
    '''
	Extracts unique OrthoGroup (OG) identifiers from the specified
    bacteria_line_file and writes them to uniq_OG file.

    Parameters:
        bacteria_line_file (str): Path to a file containing lines with
            OG identifiers.
        uniq_OG (str): Path to the output file where unique OG
            identifiers will be written.

    This function reads the lines in bacteria_line_file, extracts OG
    identifiers, and collects unique ones. It then writes these unique
    identifiers to the uniq_OG file.
    '''
    uniq_ids = set()
    with open(bacteria_line_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t') 
            if len(parts) >= 1:
                identifier = parts[0]
                uniq_ids.add(identifier)
    with open(uniq_OG, 'w') as ug:
        for identifier in uniq_ids:
            ug.write(identifier + '\n')
            
    print("Finished. the unique OG identifiers are written here ",uniq_OG)

def load_gene_names(OGs_tab_file):
    """
    Loads gene names from the provided OGs_tab_file into a dictionary.

    Parameters:
        OGs_tab_file (str): Path to the file containing OG IDs and gene names.

    Returns:
        dict: A dictionary mapping OG IDs to gene names.
    """
    gene_names = {}
    with open(OGs_tab_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                OG_id = fields[0]
                gene_name = fields[2]
                gene_names[OG_id] = gene_name
    return gene_names
    
def file_creation(uniq_OG, OrthoDB_file, final_output, OGs_tab_file):
    '''
    Creates a final output file that contains information about bacterial
    OrthoGroups (OG) based on provided data.

    Parameters:
        uniq_OG (str): Path to a file containing unique OG identifiers.
        OrthoDB_file (str): Path to the OrthoDB file.
        final_output (str): Path to the final output file.
        OGs_tab_file (str): Path to the file containing OG IDs and gene names.

    This function reads data from the uniq_OG file, the OrthoDB_file, and the
    OGs_tab_file to generate a final output file that includes OG identifiers,
    gene IDs, gene names, and species IDs. It organizes the data in a structured format.
    '''
    # Load gene names from the OGs_tab_file
    gene_names = load_gene_names(OGs_tab_file)

    with open(uniq_OG, "r") as f:
        ids = [line.strip() for line in f]

    id_to_values = {}

    with open(OrthoDB_file, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) == 2:
                id_value = fields[0]
                value = fields[1]
                if id_value in id_to_values:
                    id_to_values[id_value].append(value)
                else:
                    id_to_values[id_value] = [value]

    with open(final_output, "w") as fo:
        fo.write("OG_id\tgene_id\tspecies_id\tgene_name\n")
        for id in ids:
            values = id_to_values.get(id, [])
            if values:
                values_col2 = ';'.join(values)
                values_col3 = ';'.join(value.split(':')[0].split('_')[0] for value in values)
                gene_name = gene_names.get(id, "")
                fo.write(f"{id}\t{values_col2}\t{values_col3}\t{gene_name}\n")
            else:
                fo.write(f"{id}\t\t\t\n")
    print("Finished. The final file containing one line per bacterial OG is here :",final_output)

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################


def main():
    parser = argparse.ArgumentParser(
        description="The purpose of this script is to filter data in OrthoDB to retrieve only bacteria and format an output file with the header: OG_id gene_id species_id gene_name\
                     ",
        epilog="Exemple: python formatting_bacterial_orthologue_file.py -o ../Orthodb/odb11v0_OG2genes.tab -b only_line_bacteria.txt -u uniq_og_ids.txt -s ../Orthodb/odb11v0_level2species.tab -g ../Orthodb/odb11v0_OGs.tab -f Bacterial_OG.tab"
            )
    parser.add_argument('-o','--orthoDB_file',dest="OrthoDB_file", help="INPUT: odb11v0_OG2genes.tab",required=True)
    parser.add_argument('-g','--OGs_tab_file', dest="OGs_tab_file", help="INPUT: file containing OG IDs and gene names : odb11v0_OGs.tab",required=True)
    parser.add_argument('-s','--species_file', dest="species_file", help="INPUT: file containing bacterial species identifiers : odb11v0_level2species.tab",required=True)
    parser.add_argument('-b','--bacteria_line_file', dest="bacteria_line_file", help="OUTPUT: file with lines from odb11v0_OG2genes.tab corresponding only to the kingdom bacteria",required=True)
    parser.add_argument('-u','--uniq_OG', dest="uniq_OG", help="OUTPUT: file containing unduplicated ortholog identifiers",required=True)
    parser.add_argument('-f','--final_output', dest="final_output", help="OUTPUT: final output containing one line per bacterial OG",required=True)
    args = parser.parse_args()

    try:
        extract_line_bacteria(args.OrthoDB_file, args.bacteria_line_file, args.species_file)
        extract_unique_og_ids(args.bacteria_line_file, args.uniq_OG)
        file_creation(args.uniq_OG, args.OrthoDB_file, args.final_output, args.OGs_tab_file)
    except Exception as e:
        print(f"an error has occured : {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
