#!/usr/bin/env python3

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

def run_curl_command(protein_id):
    """
    Run the curl command to obtain information about the proteins in the OG.
    This retrieves the cds_id used to download the FASTA.
    """
    curl_command = f"curl 'https://data.orthodb.org/current/ogdetails?id={protein_id}'"
    result = subprocess.run(curl_command, shell=True, capture_output=True, text=True)
    return result

def download_fasta_content(emblcds_id):
    """
    Run the EBI API command to obtain the FASTA.
    """
    fasta_url = f'https://www.ebi.ac.uk/ena/browser/api/fasta/{emblcds_id}'
    response = requests.get(fasta_url)
    return response.text if response.status_code == 200 else None

def process_protein(protein_id):
    """
    Process a single protein ID:
      - Call OrthoDB (curl) to find its EMBLCDS reference
      - Use the EMBLCDS to fetch the FASTA from EBI
      - Build a FASTA-formatted string
    Returns (fasta_string, log_information).
    """
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
                        # 1-second pause to avoid API cancellation
                        time.sleep(1)
                        return protein_fasta, log_info
                else:
                    log_info.append(f"For ID {protein_id}, no 'EMBLCDS' ID found.")
            else:
                log_info.append(f"For ID {protein_id}, no 'xrefs' data found.")
        return None, log_info
    except requests.exceptions.SSLError as e:
        error_message = ("requests.exceptions.SSLError: None: Max retries exceeded with url: "
                         "this error comes from the API, please try again")
        return None, [error_message]

def process_row(row, processed_ogs):
    """
    Process an entire row from the original TSV:
      - Extract OG_ID
      - Retrieve the FASTAs for each protein in that OG
      - Write the FASTA and log files
      - Count the number of sequences in the resulting FASTA file
      - Return (og_id, num_sequences) or None if it fails
    """
    og_id = row['OG_ID']
    protein_ids = row['ProteinID'].split(';')

    log_filename = f'{og_id}_log.txt'
    fasta_filename = f'{og_id}_fasta.fa'

    if og_id in processed_ogs:
        return None

    # Use a temporary file for writing FASTA content.
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_fasta_file:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = list(executor.map(process_protein, protein_ids))
        for fasta_data, _ in results:
            if fasta_data:
                temp_fasta_file.write(fasta_data)

    # If the content has been written successfully, rename the temp file to the final file.
    if os.path.isfile(temp_fasta_file.name) and os.path.getsize(temp_fasta_file.name) > 0:
        shutil.move(temp_fasta_file.name, fasta_filename)
    else:
        # If the temp file is empty or wasn't created, remove it
        os.remove(temp_fasta_file.name)
        # Add OG to processed set to avoid re-processing
        processed_ogs.add(og_id)
        return None

    # Write the log file
    with open(log_filename, 'w') as log_file:
        for _, log_info in results:
            log_file.write('\n'.join(log_info) + '\n')

    # Count how many sequences in the resulting FASTA
    num_sequences = count_sequences_in_fasta(fasta_filename)
    return og_id, num_sequences

def count_sequences_in_fasta(fasta_filename):
    """
    Count the number of sequences in a FASTA file.
    """
    with open(fasta_filename, 'r') as fasta_file:
        return sum(1 for line in fasta_file if line.startswith('>'))


def fetch_orthodb_data(og_id):
    """
    Retrieves JSON information for a given OG_ID from the OrthoDB (v12) API.
    Returns a Python dictionary (or None if there's an error).
    """
    url = f"https://data.orthodb.org/v12/group?id={og_id}"
    try:
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            return resp.json().get("data", {})
        else:
            print(f"Error: Impossible to retrieve {og_id}, HTTP code = {resp.status_code}")
            return None
    except requests.exceptions.RequestException as e:
        print(f"Exception while fetching JSON for {og_id}: {e}")
        return None

def create_html_tabs(data_per_og):
    """
    Receives a dictionary data_per_og = { og_id: { 'tsv_row': {...}, 'json_data': {...} }, ... }
    and constructs an HTML string containing:
      - A list of tabs (one per OG_ID)
      - The detailed content for each tab.
    """
    html_output = """
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Genes information summary</title>
  <style>
    body {
      font-family: Arial, sans-serif;
      margin: 20px;
    }

    h1 {
      text-align: center;
    }

    .tab-button {
      background-color: #F0B64E; 
      border: 1px solid #ccc;
      padding: 10px;
      cursor: pointer;
      display: inline-block;
      margin-right: 5px;
      margin-bottom: 10px;
    }
    .tab-button:hover {
      background-color: #e5ab45; 
    }

    .tab-content {
      display: none; 
      border: 1px solid #ccc;
      padding: 10px;
      margin-bottom: 15px;
    }

    .active {
      display: block; 
    }

    table {
      border-collapse: collapse;
      margin-top: 10px;
      width: 100%;
      max-width: 800px;
      margin-bottom: 10px;
    }

    table, th, td {
      border: 1px solid #888;
      padding: 5px;
    }

    th {
      background-color: #F0B64E; 
    }

    .section-title {
      color: #B33F62;
      font-size: 1.3em;
      font-weight: bold;
      margin-top: 20px;
      margin-bottom: 10px;
    }
    .subheading {
      margin-left: 20px;
      font-weight: bold;
      margin-top: 10px;
      margin-bottom: 5px;
    }

    .centered {
      text-align: center;
    }

    .gene-name {
      color: #B33F62;
    }
  </style>
  <script>
    function openTab(evt, tabId) {
      // Hide all tab contents
      var tabcontents = document.getElementsByClassName("tab-content");
      for (var i = 0; i < tabcontents.length; i++) {
          tabcontents[i].style.display = "none";
      }
      // Show the content for the clicked tab
      document.getElementById(tabId).style.display = "block";
    }
  </script>
</head>
<body>
  <h1>Genes information summary</h1>
  <div>
"""

    # Create one tab button per OG_ID
    for i, og_id in enumerate(data_per_og.keys()):
        html_output += f"""    <button class="tab-button" onclick="openTab(event, '{og_id}')">{og_id}</button>\n"""
    html_output += "  </div>\n"

    # Build content for each tab
    for og_id, content in data_per_og.items():
        row = content['tsv_row']   # from the updated TSV
        json_data = content['json_data']  # from the OrthoDB API

        # Name of the protein (from JSON)
        protein_name = json_data.get("name", "No_name_in_JSON")

        # OG information (from the TSV)
        protein_count      = row.get("ProteinCount", "")
        species_count      = row.get("SpeciesCount", "")
        nb_single_copy     = row.get("nb_single_copy", "")
        perc_single_copy   = row.get("percent_single_copy", "")
        target_species_ct  = row.get("TargetSpecies_Count", "")
        target_species_pct = row.get("TargetSpecies_Percentage", "")
        number_of_seq      = row.get("NumberOfSeq", "")

        # Retrieve functional data
        functional_categories = json_data.get("functional_category", [])
        go_mf = json_data.get("molecular_function", [])
        go_cc = json_data.get("cellular_component", [])
        interpros = json_data.get("interpro_domains", [])
        ec_numbers = json_data.get("ECnumber", [])

        # Source link in OrthoDB
        source_url = f"https://www.orthodb.org/?level=&species=&query={og_id}"

        # Build the HTML for this OG
        html_output += f"""
  <div id="{og_id}" class="tab-content">
    <!-- OG ID and protein name centered -->
    <h2 class="centered">{og_id}</h2>
    <h3 class="centered gene-name">{protein_name}</h3>
    <p class="centered"><strong>Source:&nbsp;</strong>
       <a href="{source_url}" target="_blank">{source_url}</a>
    </p>

    <h3 class="section-title">OG information</h3>
    <table>
      <tr>
        <th>ProteinCount</th><th>SpeciesCount</th><th>nb_single_copy</th>
        <th>percent_single_copy</th><th>TargetSpecies_Count</th><th>TargetSpecies_Percentage</th>
        <th>NumberOfSeq</th>
      </tr>
      <tr>
        <td>{protein_count}</td><td>{species_count}</td><td>{nb_single_copy}</td>
        <td>{perc_single_copy}</td><td>{target_species_ct}</td><td>{target_species_pct}</td>
        <td>{number_of_seq}</td>
      </tr>
    </table>

    <h3 class="section-title">Functional Descriptions</h3>
"""

        # Functional Categories
        if functional_categories:
            html_output += "<p class='subheading'>Functional Category:</p>\n<ul>\n"
            for fc in functional_categories:
                desc = fc.get("description", "N/A")
                html_output += f"  <li>{desc}</li>\n"
            html_output += "</ul>\n"

        # GO Molecular Function
        if go_mf:
            html_output += "<p class='subheading'>GO Molecular Function:</p>\n<ul>\n"
            for mf in go_mf:
                desc   = mf.get("description", "")
                goid   = mf.get("id", "")
                count  = mf.get("count", "?")
                link   = f"https://www.ebi.ac.uk/QuickGO/GTerm?id={goid}"
                html_output += (
                    f"  <li>{count} genes with "
                    f"<a href='{link}' target='_blank'>{goid}</a>: {desc}</li>\n"
                )
            html_output += "</ul>\n"

        # GO Cellular Component
        if go_cc:
            html_output += "<p class='subheading'>GO Cellular Component:</p>\n<ul>\n"
            for cc in go_cc:
                desc   = cc.get("description", "")
                goid   = cc.get("id", "")
                count  = cc.get("count", "?")
                link   = f"https://www.ebi.ac.uk/QuickGO/GTerm?id={goid}"
                html_output += (
                    f"  <li>{count} genes with "
                    f"<a href='{link}' target='_blank'>{goid}</a>: {desc}</li>\n"
                )
            html_output += "</ul>\n"

        # InterPro
        if interpros:
            html_output += "<p class='subheading'>InterPro Domains:</p>\n<ul>\n"
            for ipr in interpros:
                desc   = ipr.get("description", "")
                ipr_id = ipr.get("id", "")
                count  = ipr.get("count", "?")
                link   = f"http://www.ebi.ac.uk/interpro/entry/InterPro/{ipr_id}/"
                html_output += (
                    f"  <li>{count} genes with "
                    f"<a href='{link}' target='_blank'>{ipr_id}</a>: {desc}</li>\n"
                )
            html_output += "</ul>\n"

        # EC numbers
        if ec_numbers:
            html_output += "<p class='subheading'>EC Numbers:</p>\n<ul>\n"
            for ec in ec_numbers:
                desc    = ec.get("description", "")
                ec_id   = ec.get("id", "")
                count   = ec.get("count", "?")
                link    = f"https://www.rhea-db.org/rhea?query=ec:{ec_id}"
                html_output += (
                    f"  <li>{count} genes with "
                    f"<a href='{link}' target='_blank'>{ec_id}</a>: {desc}</li>\n"
                )
            html_output += "</ul>\n"

        html_output += "  </div>\n"

    # End of HTML
    html_output += """
</body>
</html>
"""
    return html_output


################################################################################
# MAIN
################################################################################

def main():
    """
    Combined script:
      1) Reads an input TSV (OG_ID, ProteinID, etc.)
      2) For each OG, retrieves FASTAs using OrthoDB + EMBL API
      3) Updates the TSV with a new 'NumberOfSeq' column
      4) Generates an HTML report from the updated TSV
    """
    parser = argparse.ArgumentParser(
        description='retrieves nucleic FASTAs for each OG and generates an HTML report.',
        epilog="Example: python fastas_recovery.py OG_selected_1760.tab"
    )
    parser.add_argument('filename', help='TSV file containing OG selected information (step 2).')
    args = parser.parse_args()

    if not os.path.isfile(args.filename):
        print(f"File {args.filename} does not exist.")
        return

    # --- 1) FASTA retrieval process and updating the TSV with 'NumberOfSeq' ---

    processed_ogs = set()

    # Check which OGs have logs (to avoid re-processing)
    for filename in os.listdir('.'):
        if filename.endswith('_log.txt'):
            processed_ogs.add(filename.split('_')[0])

    # We'll build an updated version of the input TSV
    new_rows = []
    with open(args.filename, 'r') as og_file:
        og_reader = csv.DictReader(og_file, delimiter='\t')
        for row in og_reader:
            result = process_row(row, processed_ogs)
            if result is not None:
                og_id, num_sequences = result
                row['NumberOfSeq'] = num_sequences
            new_rows.append(row)

    updated_filename = f'updated_{os.path.basename(args.filename)}'
    with open(updated_filename, 'w', newline='') as output_file:
        fieldnames = og_reader.fieldnames + ['NumberOfSeq']
        writer = csv.DictWriter(output_file, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(new_rows)

    print(f"Updated TSV generated: {updated_filename}")

    # --- 2) Generate the HTML from the updated TSV ---

    # Read the updated TSV
    data_per_og = {}  # og_id -> { "tsv_row": {...}, "json_data": {...} }
    with open(updated_filename, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            og_id = row["OG_ID"]
            data_per_og[og_id] = {
                "tsv_row": row,
                "json_data": {}
            }

    # Fetch OrthoDB data for each OG_ID
    for og_id in data_per_og:
        json_data = fetch_orthodb_data(og_id)
        if json_data:
            data_per_og[og_id]["json_data"] = json_data
        else:
            data_per_og[og_id]["json_data"] = {}

    # Build the HTML content
    html_content = create_html_tabs(data_per_og)
    output_html = "rapport_OG.html"
    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"HTML file generated: {output_html}")

if __name__ == "__main__":
    main()
