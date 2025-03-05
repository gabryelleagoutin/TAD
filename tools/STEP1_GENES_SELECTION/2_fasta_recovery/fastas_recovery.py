#!/usr/bin/env python

import concurrent.futures
import subprocess
import json
import requests
import csv
import argparse
import os
import time
import tempfile
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

def run_curl_command(protein_id):
    '''Run the curl command to obtain information on the proteins contained in the OG. It will retrieve the cds_id for downloading the fasta'''
    curl_command = f"curl 'https://data.orthodb.org/current/ogdetails?id={protein_id}'"
    result = subprocess.run(curl_command, shell=True, capture_output=True, text=True)
    return result

def download_fasta_content(emblcds_id):
    ''' Run the EBI API command to obtain the fasta'''
    fasta_url = f'https://www.ebi.ac.uk/ena/browser/api/fasta/{emblcds_id}'
    response = requests.get(fasta_url)
    return response.text if response.status_code == 200 else None

def process_protein(protein_id):
    try:
        result = run_curl_command(protein_id)

        if result.returncode == 0:
            data = json.loads(result.stdout)

            log_info = []
            if "xrefs" in data["data"]:
                emblcds_info = next((xref for xref in data["data"]["xrefs"] if xref.get("type") == "EMBLCDS"), None)

                if emblcds_info:
                    emblcds_id = emblcds_info["id"]
                    log_info.append(f"For ID {protein_id}, 'EMBLCDS' ID is: {emblcds_id}")
                    fasta_data = download_fasta_content(emblcds_id)

                    if fasta_data:
                        taxid = protein_id.split(':')[0].split('_')[0]
                        species_info = fasta_data.split('\n')[0].split('|')[-1].strip()
                        protein_fasta = f">{emblcds_id}| taxid={taxid}; {species_info}\n"
                        protein_fasta += '\n'.join(fasta_data.split('\n')[1:])
                        time.sleep(1)  # Add a one-second delay between each download to avoid being cancelled by the API.
                        return protein_fasta, log_info
                else:
                    log_info.append(f"For ID {protein_id}, no 'EMBLCDS' ID found.")
            else:
                log_info.append(f"For ID {protein_id}, no 'xrefs' data found.")
        return None, log_info
    except requests.exceptions.SSLError as e:
        error_message = "requests.exceptions.SSLError: None: Max retries exceeded with url: this error comes from the API, please try again"
        return None, [error_message]

def process_row(row, processed_ogs):
    og_id = row['OG_ID']
    protein_ids = row['ProteinID'].split(';')

    log_filename = f'{og_id}_log.txt'
    fasta_filename = f'{og_id}_fasta.fa'

    if og_id in processed_ogs:
        return None

    # Create a temporary file to write content to. This means that if the program stops (e.g. if you get canceled by the API), you can re-learn where you left off. The file is only validated once it has been completely filled. Créer un fichier temporaire pour écrire le contenu
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_fasta_file:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = list(executor.map(process_protein, protein_ids))
        for fasta_data, _ in results:
            if fasta_data:
                temp_fasta_file.write(fasta_data)

    # If the content has been written successfully, rename the temporary file to the final file.
    if os.path.isfile(temp_fasta_file.name) and os.path.getsize(temp_fasta_file.name) > 0:
        shutil.move(temp_fasta_file.name, fasta_filename)
    else:
        # If the temporary file is empty or has not been created, delete the temporary file.
        os.remove(temp_fasta_file.name)

        # Add the OG to the list of processed OGs to avoid processing it again
        processed_ogs.add(og_id)

        # Exit function as no final file has been created
        return None

    with open(log_filename, 'w') as log_file:
        for _, log_info in results:
            log_file.write('\n'.join(log_info) + '\n')

    # Count the number of sequences in the FASTA file
    num_sequences = count_sequences_in_fasta(fasta_filename)

    return og_id, num_sequences

def count_sequences_in_fasta(fasta_filename):
    '''Count the number of sequences in a FASTA file.'''
    with open(fasta_filename, 'r') as fasta_file:
        return sum(1 for line in fasta_file if line.startswith('>'))

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

def main():
    parser = argparse.ArgumentParser(description='Next we need the nucleic sequences. This script retrieves the OGs selected in step 2. It takes the protein ID of each protein in the OG and then uses the OrthoDB API to retrieve the embl ID of the CDS. This xrefs is used to download the nucleic fasta with the embl API. If there is no id, it will be indicated in the log. ',
       epilog="Exemple: python fasta_recovery.py ../2_search_taxid_and_monocopy_calculation/OG_selected_1760.tab")
    parser.add_argument('filename', help='TSV file containing OG selected information in step2')
    args = parser.parse_args()

    if not os.path.isfile(args.filename):
        print(f"The file {args.filename} does not exist.")
        return

    processed_ogs = set()

    # Recover previously processed OGs
    for filename in os.listdir('.'):
        if filename.endswith('_log.txt'):
            processed_ogs.add(filename.split('_')[0])

    new_rows = []
    with open(args.filename, 'r') as og_file:
        og_reader = csv.DictReader(og_file, delimiter='\t')
        for row in og_reader:
            result = process_row(row, processed_ogs)
            if result is not None:
                og_id, num_sequences = result
                row['NumberOfSeq'] = num_sequences
            new_rows.append(row)

    output_filename = f'updated_{os.path.basename(args.filename)}'
    with open(output_filename, 'w', newline='') as output_file:
        fieldnames = og_reader.fieldnames + ['NumberOfSeq']
        writer = csv.DictWriter(output_file, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(new_rows)

if __name__ == "__main__":
    main()
