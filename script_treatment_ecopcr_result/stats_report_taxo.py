#!/usr/bin/env python

__author__ = 'Gabryelle Agoutin - INRAE & Jean Mainguy - Genoscope'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'gabryelle.agoutin@inrae.fr'
__status__ = 'prod'


import argparse
import re
import sys

# Set of ambiguous words to be cleaned up
clean_words_set = {
    'sp.',
    'Unassigned',
    'unknown',
    'metagenome',
    'unspecified',
    'unidentified',
    'unclassified',
    'incertae'
}

def filter_lines_by_keywords(cluster_file, include_keywords=None, exclude_keywords=None, clean_ambiguous=False, all_species=False, selected_species=False):
    cluster_pattern = re.compile(r"Cluster (\d+):")
    filtered_clusters = []
    rejected_clusters = []
    current_cluster_content = {}
    current_cluster = None
    all_taxonomies = set()

    def process_cluster(current_cluster, current_cluster_content):
        '''If an action is repeated several times in a function, you can use a nested function because it won't be used anywhere other than in the function and you won't have to copy and paste.'''
        # Line cleanup if the --clean_ambiguous option is enabled
        clean_passed_lines = {}
        for seq_id, content in current_cluster_content.items():
            parts = content.split("\t")
            if len(parts) != 3:
                continue
            taxonomy = parts[2]
            if clean_ambiguous and any(word in taxonomy for word in clean_words_set):
                continue
            clean_passed_lines[seq_id] = content

        if selected_species:
            valid_lines = {}
            invalid_lines = {}
            for seq_id, content in clean_passed_lines.items():
                parts = content.split("\t")
                taxonomy = parts[2]
                include_match = not include_keywords or any(keyword in taxonomy for keyword in include_keywords)
                exclude_match = exclude_keywords and any(keyword in taxonomy for keyword in exclude_keywords)
                if include_match and not exclude_match:
                    valid_lines[seq_id] = content
                elif exclude_match:
                    invalid_lines[seq_id] = content

            #  Delete invalid lines while keeping valid lines
            valid_lines = {k: v for k, v in valid_lines.items() if k not in invalid_lines}

            if valid_lines:
                filtered_clusters.append((current_cluster, valid_lines))
            else:
                rejected_clusters.append((current_cluster, current_cluster_content))
        elif all_species:
            valid_cluster = False
            for content in clean_passed_lines.values():
                parts = content.split("\t")
                taxonomy = parts[2]
                include_match = not include_keywords or any(keyword in taxonomy for keyword in include_keywords)
                exclude_match = exclude_keywords and any(keyword in taxonomy for keyword in exclude_keywords)
                if include_match and not exclude_match:
                    valid_cluster = True
                    break
            if valid_cluster:
                filtered_clusters.append((current_cluster, clean_passed_lines))
            else:
                rejected_clusters.append((current_cluster, current_cluster_content))
        else:
            filtered_lines = {}
            for seq_id, content in clean_passed_lines.items():
                parts = content.split("\t")
                taxonomy = parts[2]
                include_match = not include_keywords or any(keyword in taxonomy for keyword in include_keywords)
                exclude_match = exclude_keywords and any(keyword in taxonomy for keyword in exclude_keywords)
                if include_match and not exclude_match:
                    filtered_lines[seq_id] = content
            if filtered_lines:
                filtered_clusters.append((current_cluster, filtered_lines))
            else:
                rejected_clusters.append((current_cluster, current_cluster_content))

    with open(cluster_file, 'r') as f:
        for line in f:
            line = line.strip()
            cluster_match = cluster_pattern.match(line)
            if cluster_match:
                if current_cluster is not None:
                    process_cluster(current_cluster, current_cluster_content)
                current_cluster = cluster_match.group(1)
                current_cluster_content = {}
            else:
                parts = line.split("\t")
                if len(parts) == 3:
                    seq_id = parts[0].strip()
                    current_cluster_content[seq_id] = line

        # Process the last cluster
        if current_cluster is not None:
            process_cluster(current_cluster, current_cluster_content)

    return filtered_clusters, rejected_clusters, all_taxonomies

def calculate_statistics(filtered_clusters, output_stats_file, output_taxo_file, output_taxo_good_discriminated_file):
    total_clusters = len(filtered_clusters)
    good_discrimination = 0
    bad_discrimination = 0
    bad_clusters = []
    unique_taxonomies = set()
    good_discriminated_taxonomies = set()
    bad_discriminated_taxonomies = set()
    taxonomy_good_clusters = {}
    taxonomy_bad_clusters = {}

    for cluster, content in filtered_clusters:
        taxonomies = set()
        for line in content.values():
            parts = line.split("\t")
            if len(parts) == 3:
                taxonomy = parts[2]
                taxonomies.add(taxonomy)
                unique_taxonomies.add(taxonomy)

        if len(taxonomies) == 1:
            good_discrimination += 1
            good_discriminated_taxonomies.update(taxonomies)
            taxonomy = next(iter(taxonomies))
            taxonomy_good_clusters.setdefault(taxonomy, []).append(cluster)
        else:
            bad_discrimination += 1
            bad_clusters.append((cluster, content))
            bad_discriminated_taxonomies.update(taxonomies)
            for taxonomy in taxonomies:
                taxonomy_bad_clusters.setdefault(taxonomy, []).append(cluster)

    overlapping_taxonomies = good_discriminated_taxonomies & bad_discriminated_taxonomies
    num_overlapping_taxonomies = len(overlapping_taxonomies)

    percentage_good_discrimination = (good_discrimination / total_clusters * 100) if total_clusters > 0 else 0
    num_unique_taxonomies = len(unique_taxonomies)
    num_good_discriminated_taxonomies = len(good_discriminated_taxonomies)

    percentage_good_taxonomies = (num_good_discriminated_taxonomies / num_unique_taxonomies * 100) if num_unique_taxonomies > 0 else 0

    with open(output_stats_file, 'w') as f:
        f.write(f"Number of unique taxonomies: {num_unique_taxonomies}\n")
        f.write(f"Number of good discriminated taxonomies: {num_good_discriminated_taxonomies}\n")
        f.write(f"Percentage of good discriminated taxonomies: {percentage_good_taxonomies:.2f}%\n")
        f.write(f"Number of taxonomies both good and bad discriminated: {num_overlapping_taxonomies}\n")
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

        f.write("Taxonomies both good and bad discriminated:\n")
        for taxonomy in overlapping_taxonomies:
            f.write(f"{taxonomy}:\n")
            f.write("  Well discriminated in clusters: ")
            f.write(", ".join(taxonomy_good_clusters[taxonomy]))
            f.write("\n")
            f.write("  Poorly discriminated in clusters: ")
            f.write(", ".join(taxonomy_bad_clusters[taxonomy]))
            f.write("\n\n")

    with open(output_taxo_good_discriminated_file, 'w') as f:
        for taxonomy in sorted(good_discriminated_taxonomies):
            f.write(f"{taxonomy}\n")

    with open(output_taxo_file, 'w') as f:
        for taxonomy in sorted(unique_taxonomies):
            f.write(f"{taxonomy}\n")

    print(f"Number of unique taxonomies: {num_unique_taxonomies}")
    print(f"Number of good discriminated taxonomies: {num_good_discriminated_taxonomies}")
    print(f"Percentage of good discriminated taxonomies: {percentage_good_taxonomies:.2f}%")
    print(f"Number of taxonomies both good and bad discriminated: {num_overlapping_taxonomies}")

def generate_html_output(cluster_corrected_file, output_html_file):
    with open(cluster_corrected_file, 'r') as txt_file, open(output_html_file, 'w') as html_file:
        html_file.write("<html><body>\n")
        in_cluster = False

        for line in txt_file:
            if line.startswith("Cluster"):
                if in_cluster:
                    html_file.write("</ul>\n")
                html_file.write(f"<h2>{line.strip()}</h2>\n<ul>\n")
                in_cluster = True
            else:
                parts = line.strip().split("\t")
                if len(parts) == 3:
                    seq_id = parts[0]
                    taxid = parts[1]
                    taxonomy = parts[2]
                    html_file.write(f"<li>{seq_id}\t<a href='https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}' target='_blank'>{taxid}</a>\t{taxonomy}</li>\n")

        if in_cluster:
            html_file.write("</ul>\n")
        html_file.write("</body></html>\n")

def generate_html_output_from_stats(stats_file, output_html_file):
    with open(stats_file, 'r') as txt_file, open(output_html_file, 'w') as html_file:
        html_file.write("<html><body>\n")
        for line in txt_file:
            line = line.strip()
            if line.startswith("Number") or line.startswith("Percentage") or line.startswith("Total"):
                html_file.write(f"<p>{line}</p>\n")
            elif line.startswith("Badly discriminated clusters:"):
                html_file.write("<h2>Badly discriminated clusters:</h2>\n")
                html_file.write("<ul>\n")
            elif line.startswith("Taxonomies both good and bad discriminated:"):
                html_file.write("</ul>\n")
                html_file.write("<h2>Taxonomies both good and bad discriminated:</h2>\n")
            elif line.startswith("Cluster"):
                html_file.write(f"<h3>{line}</h3>\n<ul>\n")
            elif line and not line.startswith("k__") and not line.startswith("Cluster"):
                parts = line.split("\t")
                if len(parts) == 3:
                    seq_id = parts[0]
                    taxid = parts[1]
                    taxonomy = parts[2]
                    html_file.write(f"<li>{seq_id}\t<a href='https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}' target='_blank'>{taxid}</a>\t{taxonomy}</li>\n")
                else:
                    html_file.write(f"<li>{line}</li>\n")
            elif line.startswith("k__"):
                taxonomy = line.split(":")[0]
                html_file.write(f"<p><b>{taxonomy}:</b></p>\n<ul>\n")
            elif line.startswith("  Well discriminated"):
                html_file.write(f"<li>{line}</li>\n")
            elif line.startswith("  Poorly discriminated"):
                html_file.write(f"<li>{line}</li>\n")
            else:
                html_file.write("</ul>\n")
        html_file.write("</body></html>\n")

def main():
    parser = argparse.ArgumentParser(description="Filter clusters with optional cleaning and all_species or selected_species options.")
    parser.add_argument("-c", "--cluster_file", type=str, required=True, help="Input file containing cluster information.")
    parser.add_argument("-o", "--output_stats_file", type=str, required=True, help="Output file for the statistics report.")
    parser.add_argument("-r", "--output_corrected_file", type=str, required=True, help="Output file for the filtered clusters.")
    parser.add_argument("-i", "--include_keywords", type=str, nargs='*', help="Keywords that must be present in taxonomy to keep a line.")
    parser.add_argument("-e", "--exclude_keywords", type=str, nargs='*', help="Keywords that, if present, will cause a line to be excluded.")
    parser.add_argument("-l", "--log_file", type=str, help="Output file to log rejected clusters (optional).")
    parser.add_argument("--clean_ambiguous", action='store_true', help="Clean clusters based on a predefined dictionary of words.")
    parser.add_argument("--all_species", action='store_true', help="If set, keeps the entire cluster if it contains at least one valid line.")
    parser.add_argument("--selected_species", action='store_true', help="If set, removes clusters with only excluded lines but keeps clusters with mixed lines, removing only excluded lines.")

    args = parser.parse_args()

    if (args.include_keywords or args.exclude_keywords) and not (args.all_species or args.selected_species):
        raise ValueError("You must use --all_species or --selected_species when using -i or -e.")

    if args.all_species and args.selected_species:
        raise ValueError("You cannot use --all_species and --selected_species at the same time.")

    try:
        filtered_clusters, rejected_clusters, all_taxonomies = filter_lines_by_keywords(
            args.cluster_file,
            include_keywords=args.include_keywords,
            exclude_keywords=args.exclude_keywords,
            clean_ambiguous=args.clean_ambiguous,
            all_species=args.all_species,
            selected_species=args.selected_species
        )

        with open(args.output_corrected_file, 'w') as f:
            for cluster, content in filtered_clusters:
                f.write(f"Cluster {cluster}:\n")
                for line in content.values():
                    f.write(f"{line}\n")

        if args.log_file:
            with open(args.log_file, 'w') as log:
                for cluster, content in rejected_clusters:
                    log.write(f"Cluster {cluster}:\n")
                    for line in content.values():
                        log.write(f"{line}\n")

        calculate_statistics(
            filtered_clusters,
            args.output_stats_file,
            "uniq_taxo.txt",
            "uniq_taxo_good_discriminated.txt"
        )

        output_html_stats_file = args.output_stats_file.replace('.txt', '.html')
        generate_html_output_from_stats(args.output_stats_file, output_html_stats_file)
        output_html_file = args.output_corrected_file.replace('.txt', '.html')
        generate_html_output(args.output_corrected_file, output_html_file)

        if args.log_file:
            output_html_file_rejected = args.log_file.replace('.txt', '.html')
            generate_html_output(args.log_file, output_html_file_rejected)

        print(f"Statistics written to {args.output_stats_file}")
        print(f"Filtered clusters written to {args.output_corrected_file}")
        if args.log_file:
            print(f"Rejected clusters logged to {args.log_file}")

    except ValueError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
