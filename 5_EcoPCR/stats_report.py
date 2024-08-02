#!/usr/bin/env python

import argparse
import re

__author__ = 'Gabryelle Agoutin - INRAE'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'gabryelle.agoutin@inrae.fr'
__status__ = 'prod'

def filter_lines_by_keywords(cluster_file, taxonomy_file, output_corrected_file, include_keywords=None, exclude_keywords=None, log_file=None):
    """Filter lines within clusters based on required and excluded keywords in taxonomy."""
    
    # Compile regex patterns
    cluster_pattern = re.compile(r"Cluster (\d+):")
    taxonomy_pattern = re.compile(r"taxid=\d+; (.*)")
    
    filtered_clusters = []
    rejected_clusters = []
    current_cluster_content = {}
    current_cluster_taxonomies = set()
    current_cluster = None
    cluster_lines = {}
    
    # Parse the taxonomy file to extract taxonomy information
    taxonomy_dict = {}
    with open(taxonomy_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip()
                seq_id = header.split("|")[0][1:]
                taxonomy_dict[seq_id] = header

    # Read and process the cluster file
    with open(cluster_file, 'r') as f:
        cluster_lines = f.readlines()
        
        for line in cluster_lines:
            line = line.strip()
            cluster_match = cluster_pattern.match(line)
            if cluster_match:
                if current_cluster is not None:
                    # Filter lines within the cluster
                    filtered_lines = {}
                    for seq_id, content in current_cluster_content.items():
                        # Extract taxonomy information
                        taxonomy_match = taxonomy_pattern.search(content)
                        if taxonomy_match:
                            taxonomy = taxonomy_match.group(1)
                            if (include_keywords and any(keyword in taxonomy for keyword in include_keywords)) and \
                               (not exclude_keywords or not any(keyword in taxonomy for keyword in exclude_keywords)):
                                filtered_lines[seq_id] = content
                    # Add the cluster if it has any valid lines
                    if filtered_lines:
                        filtered_clusters.append((current_cluster, filtered_lines))
                    else:
                        rejected_clusters.append((current_cluster, current_cluster_content))
                
                # Reset for the new cluster
                current_cluster = cluster_match.group(1)
                current_cluster_taxonomies = set()
                current_cluster_content = {}
            else:
                taxonomy_match = taxonomy_pattern.search(line)
                if taxonomy_match:
                    current_cluster_taxonomies.add(taxonomy_match.group(1))
                seq_id = line.split(":")[0].strip()
                if seq_id:
                    current_cluster_content[seq_id] = line

        # Process the last cluster
        if current_cluster is not None:
            filtered_lines = {}
            for seq_id, content in current_cluster_content.items():
                taxonomy_match = taxonomy_pattern.search(content)
                if taxonomy_match:
                    taxonomy = taxonomy_match.group(1)
                    if (include_keywords and any(keyword in taxonomy for keyword in include_keywords)) and \
                       (not exclude_keywords or not any(keyword in taxonomy for keyword in exclude_keywords)):
                        filtered_lines[seq_id] = content
            if filtered_lines:
                filtered_clusters.append((current_cluster, filtered_lines))
            else:
                rejected_clusters.append((current_cluster, current_cluster_content))

    # Write the filtered clusters to a file
    with open(output_corrected_file, 'w') as f:
        for cluster, content in filtered_clusters:
            f.write(f"Cluster {cluster}:\n")
            for line in content.values():
                f.write(f"{line}\n")
            f.write("\n")

    # Write the rejected clusters if log_file is provided
    if log_file:
        with open(log_file, 'w') as log:
            log.write("Rejected clusters:\n")
            for cluster, content in rejected_clusters:
                log.write(f"Cluster {cluster}:\n")
                for line in content.values():
                    log.write(f"{line}\n")
                log.write("\n")

    return filtered_clusters

def calculate_statistics(filtered_clusters, output_stats_file):
    """Calculate and write statistics based on filtered clusters."""
    
    total_clusters = len(filtered_clusters)
    good_discrimination = 0
    bad_discrimination = 0
    bad_clusters = []

    for cluster, content in filtered_clusters:
        taxonomies = set()
        for line in content.values():
            taxonomy_match = re.search(r"taxid=\d+; (.*)", line)
            if taxonomy_match:
                taxonomies.add(taxonomy_match.group(1))
        
        if len(taxonomies) == 1:
            good_discrimination += 1
        else:
            bad_discrimination += 1
            bad_clusters.append((cluster, content))
    
    percentage_good_discrimination = (good_discrimination / total_clusters * 100) if total_clusters > 0 else 0

    # Write the statistics report
    with open(output_stats_file, 'w') as f:
        f.write(f"Total number of filtered clusters: {total_clusters}\n")
        f.write(f"Number of clusters with good discrimination: {good_discrimination}\n")
        f.write(f"Number of clusters with bad discrimination: {bad_discrimination}\n")
        f.write(f"Percentage of clusters with good discrimination: {percentage_good_discrimination:.2f}%\n\n")
        f.write("Badly discriminated clusters:\n")
        for cluster, content in bad_clusters:
            f.write(f"Cluster {cluster}:\n")
            for line in content.values():
                f.write(f"{line}\n")
            f.write("\n")

def main():
    parser = argparse.ArgumentParser(description="Filter lines within clusters based on taxonomy, generate corrected clusters and statistics.",
       epilog="python stats_report.py -c cluster.txt -t all_modified.fna -o stats.txt -i f__Lactobacillaceae -e g__Leuconostoc g__Pediococcus g__Oenococcus g__Unassigned sp. g__Weissella g__Dellaglioa g__Periweissella g__Convivina -l rejected_clusters.txt -r cluster_corrected.txt")
    parser.add_argument("-c", "--cluster_file", type=str, required=True, help="Input file containing cluster information (e.g., cluster.txt).")
    parser.add_argument("-t", "--taxonomy_file", type=str, required=True, help="Input FASTA file containing taxonomy information (e.g., all_modified.fna).")
    parser.add_argument("-o", "--output_stats_file", type=str, required=True, help="Output file for the statistics report (e.g., stats.txt).")
    parser.add_argument("-r", "--output_corrected_file", type=str, required=True, help="Output file for the filtered clusters (e.g., cluster_corrected.txt).")
    parser.add_argument("-i", "--include_keywords", type=str, nargs='*', help="Keywords that must be present in taxonomy to keep a line (e.g., f__Lactobacillaceae).")
    parser.add_argument("-e", "--exclude_keywords", type=str, nargs='*', help="Keywords that, if present, will cause a line to be excluded (e.g., f__Planococcaceae).")
    parser.add_argument("-l", "--log_file", type=str, help="Output file to log rejected clusters (optional).")
    args = parser.parse_args()

    # Filter lines within clusters based on inclusion and exclusion keywords
    filtered_clusters = filter_lines_by_keywords(
        args.cluster_file,
        args.taxonomy_file,
        args.output_corrected_file,
        include_keywords=args.include_keywords,
        exclude_keywords=args.exclude_keywords,
        log_file=args.log_file
    )

    # Generate the statistics report based on the filtered clusters
    calculate_statistics(
        filtered_clusters,
        args.output_stats_file
    )

    print(f"Statistics written to {args.output_stats_file}")
    print(f"Filtered clusters written to {args.output_corrected_file}")
    if args.log_file:
        print(f"Rejected clusters logged to {args.log_file}")

if __name__ == "__main__":
    main()
