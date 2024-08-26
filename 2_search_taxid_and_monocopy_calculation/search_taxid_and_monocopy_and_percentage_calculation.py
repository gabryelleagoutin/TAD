#!/usr/bin/env python

import argparse
import subprocess
import csv
import gzip
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import os
import shutil

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

def parse_OG_file(identifiers_with_searchID_in_taxonomy, input_file_path, taxid_to_species, min_genomes_threshold=1):
    """
    Parse the Orthologous Groups (OG) file, extracting relevant information based on specified criteria.

    :param identifiers_with_searchID_in_taxonomy: Set of identifiers with desired ID in their taxonomy
    :param input_file_path: Path to the input OG file with all bacterian OG
    :param taxid_to_species: Dictionary mapping taxonomy IDs to species names.
    :param min_genomes_threshold: Minimal number of genomes for selecting OGs (default is 1).
    :return: List of dictionaries containing extracted information for each OG.

    The OG file is tab-delimited and expected to have the following columns:
    - Column 1: OG_ID (Orthologous Group ID)
    - Column 2: ProteinID
    - Column 3: SpeciesID (semicolon-separated list of taxonomy IDs)

    For each OG that matches the specified identifiers (if a taxid contained in identifiers_with_searchID_in_taxonomy is present in the speciesID column of the OG, we take the whole OG)
     the function calculates and includes the following information:
    - 'OG_ID': Orthologous Group ID
    - 'ProteinCount': Total number of proteins in the OG
    - 'SpeciesCount': Number of unique species in the OG
    - 'nb_single_copy': Number of proteins present as a single copy in at least one species
    - 'percent_single_copy': Percentage of proteins that are single-copy in the OG
    - 'ProteinID': ID of a representative protein in the OG
    - 'taxids': List of taxonomy IDs corresponding to the species in the OG
    - 'species': List of species names corresponding to the taxonomy IDs
    - 'TargetSpecies_Count': Number of target species in the OG
    - 'TargetSpecies_Percentage': Percentage of target species in the OG

    The function returns a list of dictionaries, each containing information for a selected OG.
    If no OGs meet the criteria, a ValueError is raised.
    """

    cogs = []
    proper_open = gzip.open if input_file_path.endswith('.gz') else open
    taxids_set = set()

    with proper_open(input_file_path, 'r') as input_file:
        for line in input_file:
            columns = line.strip().split('\t')
            if len(columns) >= 3:
                ids_in_column = columns[2].split(';')
                matching_ids = [id for id in identifiers_with_searchID_in_taxonomy if id in ids_in_column]
                if matching_ids:
                    OG_ID, ProteinID, SpeciesID ,GeneName = line.strip().split('\t')
                    tax_ids = [taxid.split(':')[0] for taxid in ProteinID.split(';')]
                    sp_count = len(tax_ids)

                    if sp_count >= min_genomes_threshold:
                        prot_count = 0
                        copy_counts = {}
                        for protid in tax_ids:
                            if protid in copy_counts:
                                copy_counts[protid] += 1
                            else:
                                copy_counts[protid] = 1

                        new_taxid = set(SpeciesID.split(';'))

                        nb_single_copy = sum((1 for gene_copy in copy_counts.values() if gene_copy == 1))

                        target_species_count = sum(1 for taxid in new_taxid if taxid in identifiers_with_searchID_in_taxonomy)
                        target_species_percentage = (target_species_count / len(new_taxid)) * 100 if new_taxid else 0

                        cog = {'OG_ID': OG_ID,
                               'ProteinCount': int(sp_count),
                               'SpeciesCount': int(len(new_taxid)),
                               'nb_single_copy': int(nb_single_copy),
                               'percent_single_copy': (int(nb_single_copy) / int(len(copy_counts))) * 100,
                               'ProteinID': ProteinID,
                               'taxids': new_taxid,
                               'species': map_species(new_taxid, taxid_to_species),
                               'gene_name': GeneName,
                               'TargetSpecies_Count': target_species_count,
                               'TargetSpecies_Percentage': round(target_species_percentage, 2)
                               }

                        cogs.append(cog)

    if len(cogs) == 0:
        raise ValueError('No COGs identified due to the min_genome_threshold')
    return cogs

def filter_matching_lines(input_file_path, search_ID, level2species_path):
    """
    Search for the IDs of species that have the desired ID in their taxonomy.
    """
    command = f"grep -w {search_ID} {level2species_path} | cut -f2"
    result = subprocess.check_output(command, shell=True, text=True)

    identifiers_with_searchID_in_taxonomy = set(result.strip().split('\n'))
    identifiers_with_searchID_in_taxonomy = {id.split('_')[0] for id in identifiers_with_searchID_in_taxonomy}
    return identifiers_with_searchID_in_taxonomy

def map_species(taxids, taxid_to_species):
    """
    Lists the NAMES of the species in the table.
    """
    species_list = [taxid_to_species.get(taxid, '') for taxid in taxids]
    return species_list

##################################################################################################################################################
#
# PARALLELISATION
#
##################################################################################################################################################

def process_file(input_file, search_ID, taxid_to_species, min_genomes_threshold, output_file, level2species_path):
    """
    Processes the input file in chunks, searching for identifiers with the specified ID in their taxonomy and calling process_input_file.
    """
    identifiers_with_searchID_in_taxonomy = filter_matching_lines(input_file, search_ID, level2species_path)
    process_input_file({'input_file': input_file, 'search_ID': search_ID, 'min_genomes_threshold': min_genomes_threshold, 'taxid_to_species': taxid_to_species, 'output_file': output_file}, identifiers_with_searchID_in_taxonomy)

def process_input_file(args, identifiers_with_searchID_in_taxonomy):
    """
    Processes the input file in chunks, calling the parse_eggnog_members_file function.
    """
    taxid_to_species = args['taxid_to_species']
    input_file = args['input_file']
    search_ID = args['search_ID']
    output_file = args['output_file']

    if not identifiers_with_searchID_in_taxonomy or (len(identifiers_with_searchID_in_taxonomy) == 1 and '' in identifiers_with_searchID_in_taxonomy):
        raise ValueError(f'The number {search_ID} was not found in the file.')

    result = parse_OG_file(identifiers_with_searchID_in_taxonomy, input_file, taxid_to_species, args['min_genomes_threshold'])

    with open(output_file, 'a+', newline='') as tsvfile:
        fieldnames = ['OG_ID', 'ProteinCount', 'SpeciesCount', 'nb_single_copy', 'percent_single_copy', 'ProteinID', 'taxids', 'species', 'gene_name', 'TargetSpecies_Count', 'TargetSpecies_Percentage']
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')

        for cog in result:
            cog['taxids'] = ', '.join([taxid.strip('"') for taxid in cog['taxids']])
            writer.writerow(cog)

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

def main():
    parser = argparse.ArgumentParser(
        description="""For each OG that matches the specified identifiers the script calculates and includes:\
    - 'OG_ID': Orthologous Group ID\
    - 'ProteinCount': Total number of proteins in the OG\
    - 'SpeciesCount': Number of unique species in the OG\
    - 'nb_single_copy': Number of proteins present as a single copy in at least one species\
    - 'percent_single_copy': Percentage of proteins that are single-copy in the OG\
    - 'ProteinID': List of all protein IDs. Protein IDs are made up of a number followed by ‘:’ followed by the taxid.\
    - 'taxids': List of taxonomy IDs corresponding to the species in the OG\
    - 'species': List of species names corresponding to the taxonomy IDs\
    - 'TargetSpecies_Count': Number of target species in the OG\
    - 'TargetSpecies_Percentage': Percentage of target species in the OG""",
       epilog="Exemple: python search_taxid_and_monocopy_calculation.py -i Bacterial_OG.tab -f ../Orthodb/odb11v0_species.tab -s 1578 -l ../Orthodb/odb11v0_level2species.tab -o OG_1578.tab")
    parser.add_argument("-i", "--input_file", help="file containing all bacterial orthologue groups.is tab-delimited and expected to have the following columns: OG_ID,ProteinID,speciesID")
    parser.add_argument("-s", "--search_ID", type=int, help="Identifier of the taxonomic rank you are looking for. Example: for Lactobacillus, the identifier is 1568")
    parser.add_argument('--min_genomes_threshold', type=int, default=1, help='Minimal number of genomes for cog selection')
    parser.add_argument("-o", '--output_tsv', default='OG_stat_single_copy.tsv', help='Path to the output TSV file')
    parser.add_argument("-f", '--species_file', help='Path to the species file (e.g., odb11v0_species.tab)')
    parser.add_argument("-l", '--level2species_file', help='Path to the level2species file (e.g., odb11v0_level2species.tab)', required=True)
    args = parser.parse_args()

    if not args.species_file:
        parser.error("The species file (-f/--species_file) is required.")

    with open(args.species_file, 'r') as file:
        odb11v0_species = pd.read_csv(file, delimiter='\t', header=None, names=['NCBI_taxid', 'orthoDB_taxid', 'species', 'genome_id', 'genome_size', 'OG_count', 'coding'])
    taxid_to_species = dict(zip(odb11v0_species['orthoDB_taxid'].str.split('_').str[0], odb11v0_species['species']))

    # Check if the temp_chunks directory exists
    temp_dir = 'temp_chunks'
    if os.path.exists(temp_dir):
        print(f"Directory '{temp_dir}' already exists. Removing existing directory...")
        # Remove the existing directory and its contents
        shutil.rmtree(temp_dir)
        print(f"Directory '{temp_dir}' removed.")

    # Create a temporary directory for storing chunks
    temp_dir = 'temp_chunks'
    os.makedirs(temp_dir, exist_ok=True)

    # Split the input file into chunks
    chunk_size = 1000
    with open(args.input_file, 'r') as input_file:
        chunk_count = 0
        chunk = []
        for line in input_file:
            if len(chunk) >= chunk_size:
                with open(os.path.join(temp_dir, f'chunk_{chunk_count}.tsv'), 'w') as chunk_file:
                    chunk_file.writelines(chunk)
                chunk_count += 1
                chunk = []
            chunk.append(line)

        # Write the remaining lines to the last chunk
        with open(os.path.join(temp_dir, f'chunk_{chunk_count}.tsv'), 'w') as chunk_file:
            chunk_file.writelines(chunk)
        chunk_count += 1

    # Process chunks in parallel
    with ProcessPoolExecutor() as executor:
        futures = []
        for i in range(chunk_count):
            input_chunk_file = os.path.join(temp_dir, f'chunk_{i}.tsv')
            output_chunk_file = os.path.join(temp_dir, f'chunk_{i}_output.tsv')
            future = executor.submit(process_file, input_chunk_file, args.search_ID, taxid_to_species, args.min_genomes_threshold, output_chunk_file, args.level2species_file)
            futures.append(future)

        # Wait for all tasks to complete
        for future in futures:
            future.result()

    # Combine the output chunks into the final result
    with open(args.output_tsv, 'w', newline='') as final_output:
        fieldnames = ['OG_ID', 'ProteinCount', 'SpeciesCount', 'nb_single_copy', 'percent_single_copy', 'ProteinID', 'taxids', 'species', 'gene_name', 'TargetSpecies_Count', 'TargetSpecies_Percentage']
        writer = csv.DictWriter(final_output, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for i in range(chunk_count):
            output_chunk_file = os.path.join(temp_dir, f'chunk_{i}_output.tsv')
            with open(output_chunk_file, 'r') as chunk_file:
                for line in chunk_file:
                    # Parse the line into a dictionary
                    cog_dict = dict(zip(fieldnames, line.strip().split('\t')))
                    # Write the dictionary to the CSV file
                    writer.writerow(cog_dict)

    # Clean up temporary directory
    for i in range(chunk_count):
        input_chunk_file = os.path.join(temp_dir, f'chunk_{i}.tsv')
        output_chunk_file = os.path.join(temp_dir, f'chunk_{i}_output.tsv')
        os.remove(input_chunk_file)
        os.remove(output_chunk_file)

    os.rmdir(temp_dir)

if __name__ == "__main__":
    main()
