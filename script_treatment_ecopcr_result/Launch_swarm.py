#!/usr/bin/env python

import subprocess
import argparse


__author__ = 'Gabryelle Agoutin - INRAE'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'gabryelle.agoutin@inrae.fr'
__status__ = 'prod'


def parse_fasta_with_taxonomy(fasta_file):
    """Parse the FASTA file to extract taxonomy information."""
    taxonomy_dict = {}
    with open(fasta_file, 'r') as f:
        seq_id = None
        for line in f:
            if line.startswith(">"):
                header = line.strip()
                seq_id = header.split("|")[0][1:]
                taxonomy_dict[seq_id] = header
    return taxonomy_dict

def run_swarm(fasta_file, output_file, threads=4, abundance=1, distance=1):
    """Run Swarm to cluster sequences."""
    cmd = [
        "swarm",
        "-t", str(threads),
        "-a", str(abundance),
        "-d", str(distance),
        "-f", fasta_file,
        "-o", output_file
    ]
    subprocess.run(cmd, check=True)

def process_swarm_output(swarm_file, taxonomy_dict, output_info_file):
    """Process the Swarm output to generate a detailed report."""
    with open(swarm_file, 'r') as clusters_file, open(output_info_file, 'w') as out_file:
        clusters = clusters_file.readlines()
        for i, cluster in enumerate(clusters, start=1):
            cluster_ids = cluster.strip().split()
            out_file.write(f"Cluster {i}:\n")
            for cluster_id in cluster_ids:
                seq_id = cluster_id.split("|")[0]
                if seq_id in taxonomy_dict:
                    info = (taxonomy_dict[seq_id]).split(" ")
                    tax = info[2]
                    taxid = info[1].split(';')[0].split('=')[1]
                    out_file.write(f"{seq_id}\t{taxid}\t{tax}\n")
                else:
                    out_file.write(f"No info found for {seq_id} in the file.\n")

def process_swarm_output_html(swarm_file, taxonomy_dict, output_info_file_html):
    """Process the Swarm output to generate a detailed HTML report."""
    with open(swarm_file, 'r') as clusters_file, open(output_info_file_html, 'w') as out_file:
        clusters = clusters_file.readlines()
        out_file.write("<html><body>\n")
        for i, cluster in enumerate(clusters, start=1):
            cluster_ids = cluster.strip().split()
            out_file.write(f"<h2>Cluster {i}:</h2>\n<ul>\n")
            for cluster_id in cluster_ids:
                seq_id = cluster_id.split("|")[0]
                if seq_id in taxonomy_dict:
                    info = (taxonomy_dict[seq_id]).split(" ")
                    tax = info[2]
                    taxid = info[1].split(';')[0].split('=')[1]
                    # Create a clickable link for the taxid
                    out_file.write(f"<li>{seq_id}\t<a href='https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}' target='_blank'>{taxid}</a>\t{tax}</li>\n")
                else:
                    out_file.write(f"<li>No info found for {seq_id} in the file.</li>\n")
            out_file.write("</ul>\n")
        out_file.write("</body></html>\n")

def generate_taxonomy_table(taxonomy_dict, output_table_file="taxonomy_rank_table.txt"):
    """Generate a table with taxonomy ranks (Kingdom, Phylum, etc.) and unique counts."""
    taxonomy_levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    rank_dict = {level: set() for level in taxonomy_levels}

    with open(output_table_file, 'w') as out_file:
        out_file.write("\t".join(taxonomy_levels) + "\n")

        for tax_info in taxonomy_dict.values():
            tax_ranks = tax_info.split(" ")[2].split(";")
            if len(tax_ranks) == 7:  
                for i, level in enumerate(taxonomy_levels):
                    rank_dict[level].add(tax_ranks[i])

        # Write the unique counts for each rank
        unique_counts = [str(len(rank_dict[level])) for level in taxonomy_levels]
        out_file.write("\t".join(unique_counts) + "\n")

        # Write the unique names for each rank
        max_rows = max(len(rank_dict[level]) for level in taxonomy_levels)
        rank_lists = [list(rank_dict[level]) for level in taxonomy_levels]
        for i in range(max_rows):
            row = [
                rank_lists[j][i] if i < len(rank_lists[j]) else "" for j in range(len(taxonomy_levels))
            ]
            out_file.write("\t".join(row) + "\n")

    print(f"Taxonomy rank table written to {output_table_file}")

def main():
    parser = argparse.ArgumentParser(description="Pipeline to process EcoPCR output, run Swarm, and analyze results.",
       epilog="Exemple: python Launch_swarm.py -f all_modified.fna -s fichier_swarm.txt -o cluster.txt  -t 4 -a 1 -d 1")
    parser.add_argument("-f", "--fasta_file", type=str, required=True, help="The input FASTA file (all_modified.fna).")
    parser.add_argument("-s", "--swarm_output_file", type=str, required=True, help="The output file for Swarm (fichier_swarm.txt).")
    parser.add_argument("-o", "--output_info_file", type=str, required=True, help="The output info file (cluster.txt).")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use for Swarm.")
    parser.add_argument("-a", "--abundance", type=int, default=1, help="Minimum abundance for Swarm.")
    parser.add_argument("-d", "--distance", type=int, default=1, help="Distance threshold for Swarm.")
    parser.add_argument("-oh", "--output_info_file_html", type=str, required=False, help="The output HTML info file (cluster.html).")
    args = parser.parse_args()

    # Step 1: Parse and extract taxonomy information from the modified FASTA file
    print("Parsing FASTA file for taxonomy information...")
    taxonomy_dict = parse_fasta_with_taxonomy(args.fasta_file)

    # Step 2: Run Swarm
    print("Running Swarm...")
    run_swarm(args.fasta_file, args.swarm_output_file, args.threads, args.abundance, args.distance)

    # Step 3: Process Swarm output to generate the desired report
    print("Processing Swarm output...")
    process_swarm_output(args.swarm_output_file, taxonomy_dict, args.output_info_file)
    print(f"Results written to {args.output_info_file}")
    
    # Step 4: Generate taxonomy rank table
    generate_taxonomy_table(taxonomy_dict)
    
    
    # Step 4: Generate HTML file
    output_html_file = args.output_info_file.replace('.txt', '.html')
    process_swarm_output_html(args.swarm_output_file, taxonomy_dict,output_html_file)


if __name__ == "__main__":
    main()
