#!/usr/bin/env python

import subprocess
import argparse
import logging
import re
from collections import defaultdict

__author__ = 'Gabryelle Agoutin - INRAE'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'gabryelle.agoutin@inrae.fr'
__status__ = 'prod'

def normalize_id(seq_id):
    """
    Normalize sequence ID by removing a trailing '_X', where X is one or more digits.
    """
    match = re.match(r"^(.*?)(_\d+)?$", seq_id)
    if match:
        return match.group(1)  # Return the part before the optional '_X'
    return seq_id

def parse_name_seq_with_taxo(taxonomy_file):
    """Parse the name_seq_with_taxo.txt file to extract taxonomy information."""
    taxonomy_dict = {}
    with open(taxonomy_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip()
                parts = header[1:].split(" ", 1)
                seq_id = normalize_id(parts[0].split("|")[0])  # Normalize ID
                taxonomy_info = parts[1] if len(parts) > 1 else ""
                taxonomy_dict[seq_id] = taxonomy_info
    return taxonomy_dict

def parse_vsearch_clusters(vsearch_cluster_file):
    """Parse the output_vsearch_cluster.txt file to create a mapping of centroids to their sequences."""
    cluster_dict = {}
    with open(vsearch_cluster_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # Ignore empty lines
            parts = line.split('\t')
            if len(parts) == 2:  # Ensure there are two parts (centroid and sequences)
                centroid = normalize_id(parts[0].split("_")[0])  # Normalize centroid ID
                sequences = parts[1].split(',') if parts[1] else []
                cluster_dict[centroid] = [normalize_id(seq) for seq in sequences]  # Normalize sequences
            elif len(parts) == 1:  # Handle cases where there are no associated sequences
                centroid = normalize_id(parts[0])  # Normalize centroid ID
                cluster_dict[centroid] = []  # No associated sequences
    return cluster_dict

def augment_swarm_clusters(swarm_file, vsearch_cluster_file, output_swarm_file):
    """Augment Swarm clusters with sequences from vsearch clusters."""
    logging.info("Augmenting Swarm clusters with VSEARCH clusters...")
    vsearch_clusters = parse_vsearch_clusters(vsearch_cluster_file)
    
    with open(swarm_file, 'r') as f_in, open(output_swarm_file, 'w') as f_out:
        for line in f_in:
            cluster_sequences = line.strip().split()
            augmented_cluster = set(cluster_sequences)  # Use a set to avoid duplicates
            
            # For each sequence in the cluster, check if it's a centroid in the VSEARCH output
            for seq_full in cluster_sequences:
                seq_id = normalize_id(seq_full.split("|")[0])  # Normalize sequence ID
                
                # Check if this sequence has associated sequences in the VSEARCH clusters
                if seq_id in vsearch_clusters:
                    associated_sequences = vsearch_clusters[seq_id]
                    
                    # Add associated sequences to the cluster if not already present
                    for associated_seq in associated_sequences:
                        augmented_cluster.add(associated_seq + "|_1")
            
            # Write the augmented cluster to the output file
            f_out.write(' '.join(sorted(augmented_cluster)) + '\n')
    
    logging.info(f"Augmented Swarm clusters written to {output_swarm_file}")

def run_swarm(fasta_file, output_file, threads=4, abundance=1, distance=1):
    """Run Swarm to cluster sequences."""
    logging.info("Running Swarm...")
    cmd = [
        "swarm",
        "-t", str(threads),
        "-a", str(abundance),
        "-d", str(distance),
        "-f", fasta_file,
        "-o", output_file
    ]
    subprocess.run(cmd, check=True)
    logging.info(f"Swarm output written to {output_file}")

def process_swarm_output(swarm_file, taxonomy_dict, output_info_file):
    """Process the Swarm output to generate a detailed report."""
    logging.info("Processing Swarm output...")
    with open(swarm_file, 'r') as clusters_file, open(output_info_file, 'w') as out_file:
        clusters = clusters_file.readlines()
        for i, cluster in enumerate(clusters, start=1):
            cluster_ids = cluster.strip().split()
            out_file.write(f"Cluster {i}:\n")
            for cluster_id in cluster_ids:
                seq_id = normalize_id(cluster_id.split("|")[0])  # Normalize sequence ID
                if seq_id in taxonomy_dict:
                    info = taxonomy_dict[seq_id]
                    parts = info.split(' ', 1)
                    taxid_part = parts[0] if parts else ""
                    tax_part = parts[1] if len(parts) > 1 else ""
                    taxid = taxid_part.split('=')[1].split(';')[0] if 'taxid=' in taxid_part else "Unknown"
                    tax = tax_part.strip() if tax_part else "Unknown"
                    out_file.write(f"{seq_id}\t{taxid}\t{tax}\n")
                else:
                    out_file.write(f"No info found for {seq_id} in the file.\n")
    logging.info(f"Swarm processing results written to {output_info_file}")

def process_swarm_output_html(swarm_file, taxonomy_dict, output_info_file_html):
    """Process the Swarm output to generate a detailed HTML report."""
    logging.info("Generating HTML report...")
    with open(swarm_file, 'r') as clusters_file, open(output_info_file_html, 'w') as out_file:
        clusters = clusters_file.readlines()
        out_file.write("<html><body>\n")
        for i, cluster in enumerate(clusters, start=1):
            cluster_ids = cluster.strip().split()
            out_file.write(f"<h2>Cluster {i}:</h2>\n<ul>\n")
            for cluster_id in cluster_ids:
                seq_id = normalize_id(cluster_id.split("|")[0])  # Normalize sequence ID
                if seq_id in taxonomy_dict:
                    info = taxonomy_dict[seq_id]
                    parts = info.split(' ', 1)
                    taxid_part = parts[0] if parts else ""
                    tax_part = parts[1] if len(parts) > 1 else ""
                    taxid = taxid_part.split('=')[1].split(';')[0] if 'taxid=' in taxid_part else "Unknown"
                    tax = tax_part.strip() if tax_part else "Unknown"
                    # Create a clickable link for the taxid
                    out_file.write(f"<li>{seq_id}\t<a href='https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}' target='_blank'>{taxid}</a>\t{tax}</li>\n")
                else:
                    out_file.write(f"<li>No info found for {seq_id} in the file.</li>\n")
            out_file.write("</ul>\n")
        out_file.write("</body></html>\n")
    logging.info(f"HTML report written to {output_info_file_html}")

def generate_taxonomy_table(taxonomy_dict, output_table_file="taxonomy_rank_table.txt"):
    """Generate a table with taxonomy ranks (Kingdom, Phylum, etc.) and unique counts."""
    logging.info("Generating taxonomy rank table...")
    taxonomy_levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    rank_dict = {level: set() for level in taxonomy_levels}

    with open(output_table_file, 'w') as out_file:
        out_file.write("\t".join(taxonomy_levels) + "\n")

        for tax_info in taxonomy_dict.values():
            tax_ranks = tax_info.split(" ")[-1].split(";")
            if len(tax_ranks) == 7:
                for i, level in enumerate(taxonomy_levels):
                    rank_dict[level].add(tax_ranks[i].strip())

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
    logging.info(f"Taxonomy rank table written to {output_table_file}")

def main():
    parser = argparse.ArgumentParser(description="Pipeline to process EcoPCR output, augment Swarm clusters, and analyze results.",
       epilog="Exemple: python Launch_swarm.py -f all_modified.fna -s fichier_swarm.txt -o cluster.txt -t 4 -a 1 -d 1 -vsearch output_vsearch_cluster.txt -taxo name_seq_with_taxo.txt")
    parser.add_argument("-f", "--fasta_file", type=str, required=True, help="The input FASTA file (all_modified.fna).")
    parser.add_argument("-s", "--swarm_output_file", type=str, required=True, help="The output file for Swarm (fichier_swarm.txt).")
    parser.add_argument("-o", "--output_info_file", type=str, required=True, help="The output info file (cluster.txt).")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use for Swarm.")
    parser.add_argument("-a", "--abundance", type=int, default=1, help="Minimum abundance for Swarm.")
    parser.add_argument("-d", "--distance", type=int, default=1, help="Distance threshold for Swarm.")
    parser.add_argument("-oh", "--output_info_file_html", type=str, required=False, help="The output HTML info file (cluster.html).")
    parser.add_argument("-vsearch", "--vsearch_cluster_file", type=str, required=True, help="The VSEARCH cluster file (output_vsearch_cluster.txt).")
    parser.add_argument("-taxo", "--taxonomy_file", type=str, required=True, help="The taxonomy file (name_seq_with_taxo.txt).")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(message)s')

    # Step 1: Parse and extract taxonomy information from the taxonomy file
    logging.info("Parsing taxonomy file for taxonomy information...")
    taxonomy_dict = parse_name_seq_with_taxo(args.taxonomy_file)

    # Step 2: Run Swarm
    run_swarm(args.fasta_file, args.swarm_output_file, args.threads, args.abundance, args.distance)

    # Step 3: Augment Swarm clusters using VSEARCH clusters
    augmented_swarm_file = args.swarm_output_file.replace('.txt', '_complete.txt')
    augment_swarm_clusters(args.swarm_output_file, args.vsearch_cluster_file, augmented_swarm_file)

    # Step 4: Process augmented Swarm output to generate the desired report
    process_swarm_output(augmented_swarm_file, taxonomy_dict, args.output_info_file)
    logging.info(f"Results written to {args.output_info_file}")

    # Step 5: Generate taxonomy rank table
    generate_taxonomy_table(taxonomy_dict)

    # Step 6: Generate HTML file
    if args.output_info_file_html:
        process_swarm_output_html(augmented_swarm_file, taxonomy_dict, args.output_info_file_html)

if __name__ == "__main__":
    main()
