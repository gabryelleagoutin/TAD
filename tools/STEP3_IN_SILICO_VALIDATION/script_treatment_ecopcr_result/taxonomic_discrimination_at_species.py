import argparse
import re
import sys

# Dictionary of words to clean
clean_words_dict = {
    'sp.': True,
    'Unassigned': True,
    'unknown': True,
}

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################


def filter_lines_by_keywords(cluster_file, include_keywords=None, exclude_keywords=None, clean_words=False, all_species=False, selected_species=False):
    """
    Reads a cluster file and filters lines based on various conditions:
      - If clean_words is True, lines containing certain dictionary words are removed.
      - If selected_species is True, lines are filtered so that only lines containing
        include_keywords are kept (and exclude_keywords are removed).
      - If all_species is True, the entire cluster is kept if it contains at least one valid line.
    Returns:
      - A list of (cluster_number, {seq_id: line}) for the filtered clusters
      - A list of (cluster_number, {seq_id: line}) for the rejected clusters
      - A set of all taxonomy strings encountered
    """
    cluster_pattern = re.compile(r"Cluster (\d+):")
    filtered_clusters = []
    rejected_clusters = []
    current_cluster_content = {}
    current_cluster = None
    all_taxonomies = set()

    with open(cluster_file, 'r') as f:
        cluster_lines = f.readlines()

        for line in cluster_lines:
            line = line.strip()
            cluster_match = cluster_pattern.match(line)
            if cluster_match:
                # If we reach a new cluster, we process the previous cluster
                if current_cluster is not None:
                    # Clean lines if clean_words is activated
                    clean_passed_lines = {
                        seq_id: content for seq_id, content in current_cluster_content.items()
                        if not clean_words or not any(word in content.split("\t")[2] for word in clean_words_dict)
                    }

                    # Apply filtering for --selected_species
                    if selected_species:
                        valid_lines = {
                            seq_id: content for seq_id, content in clean_passed_lines.items()
                            if (not include_keywords or any(keyword in content.split("\t")[2] for keyword in include_keywords)) and
                               (not exclude_keywords or not any(keyword in content.split("\t")[2] for keyword in exclude_keywords))
                        }
                        invalid_lines = {
                            seq_id: content for seq_id, content in clean_passed_lines.items()
                            if exclude_keywords and any(keyword in content.split("\t")[2] for keyword in exclude_keywords)
                        }

                        # Remove invalid lines (with -e) while keeping valid ones
                        valid_lines = {k: v for k, v in valid_lines.items() if k not in invalid_lines}

                        if valid_lines:
                            filtered_clusters.append((current_cluster, valid_lines))
                        else:
                            rejected_clusters.append((current_cluster, current_cluster_content))

                    elif all_species:
                        valid_cluster = any(
                            (not include_keywords or any(keyword in content.split("\t")[2] for keyword in include_keywords)) and
                            (not exclude_keywords or not any(keyword in content.split("\t")[2] for keyword in exclude_keywords))
                            for content in clean_passed_lines.values()
                        )
                        if valid_cluster:
                            filtered_clusters.append((current_cluster, clean_passed_lines))
                        else:
                            rejected_clusters.append((current_cluster, current_cluster_content))
                    else:
                        # Standard line-by-line filtering
                        filtered_lines = {
                            seq_id: content for seq_id, content in clean_passed_lines.items()
                            if (not include_keywords or any(keyword in content.split("\t")[2] for keyword in include_keywords)) and
                               (not exclude_keywords or not any(keyword in content.split("\t")[2] for keyword in exclude_keywords))
                        }
                        if filtered_lines:
                            filtered_clusters.append((current_cluster, filtered_lines))
                        else:
                            rejected_clusters.append((current_cluster, current_cluster_content))

                current_cluster = cluster_match.group(1)
                current_cluster_content = {}
            else:
                parts = line.split("\t")
                if len(parts) == 3:
                    seq_id = parts[0].strip()
                    current_cluster_content[seq_id] = line

        # Process the last cluster in the file
        if current_cluster is not None:
            clean_passed_lines = {
                seq_id: content for seq_id, content in current_cluster_content.items()
                if not clean_words or not any(word in content.split("\t")[2] for word in clean_words_dict)
            }

            if selected_species:
                valid_lines = {
                    seq_id: content for seq_id, content in clean_passed_lines.items()
                    if (not include_keywords or any(keyword in content.split("\t")[2] for keyword in include_keywords)) and
                       (not exclude_keywords or not any(keyword in content.split("\t")[2] for keyword in exclude_keywords))
                }
                invalid_lines = {
                    seq_id: content for seq_id, content in clean_passed_lines.items()
                    if exclude_keywords and any(keyword in content.split("\t")[2] for keyword in exclude_keywords)
                }

                # Remove invalid lines (with -e) while keeping valid lines
                valid_lines = {k: v for k, v in valid_lines.items() if k not in invalid_lines}

                if valid_lines:
                    filtered_clusters.append((current_cluster, valid_lines))
                else:
                    rejected_clusters.append((current_cluster, current_cluster_content))
            elif all_species:
                valid_cluster = any(
                    (not include_keywords or any(keyword in content.split("\t")[2] for keyword in include_keywords)) and
                    (not exclude_keywords or not any(keyword in content.split("\t")[2] for keyword in exclude_keywords))
                    for content in clean_passed_lines.values()
                )
                if valid_cluster:
                    filtered_clusters.append((current_cluster, clean_passed_lines))
                else:
                    rejected_clusters.append((current_cluster, current_cluster_content))
            else:
                filtered_lines = {
                    seq_id: content for seq_id, content in clean_passed_lines.items()
                    if (not include_keywords or any(keyword in content.split("\t")[2] for keyword in include_keywords)) and
                       (not exclude_keywords or not any(keyword in content.split("\t")[2] for keyword in exclude_keywords))
                }
                if filtered_lines:
                    filtered_clusters.append((current_cluster, filtered_lines))
                else:
                    rejected_clusters.append((current_cluster, current_cluster_content))

    return filtered_clusters, rejected_clusters, all_taxonomies

def calculate_statistics(filtered_clusters, output_stats_file, output_taxo_file, output_taxo_good_discriminated_file):
    """
    Calculates various statistics about the filtered clusters:
      - Number of unique taxonomies
      - Number of well-discriminated taxonomies (taxonomies appearing in clusters containing only one taxonomy)
      - Number of overlapping taxonomies (taxonomies that appear both in a well-discriminated cluster and in a poorly discriminated cluster)
      - Writes these statistics to output_stats_file and prints them.
      - Writes the list of all unique taxonomies to output_taxo_file
      - Writes the list of well-discriminated taxonomies to output_taxo_good_discriminated_file
    Also returns:
      - The content needed to display "badly discriminated clusters"
      - The dictionary of good and bad clusters per taxonomy (for further usage if needed)
    """
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
            for taxonomy in taxonomies:
                if taxonomy not in taxonomy_good_clusters:
                    taxonomy_good_clusters[taxonomy] = []
                taxonomy_good_clusters[taxonomy].append(cluster)
        else:
            bad_discrimination += 1
            bad_clusters.append((cluster, content))
            bad_discriminated_taxonomies.update(taxonomies)
            for taxonomy in taxonomies:
                if taxonomy not in taxonomy_bad_clusters:
                    taxonomy_bad_clusters[taxonomy] = []
                taxonomy_bad_clusters[taxonomy].append(cluster)

    overlapping_taxonomies = good_discriminated_taxonomies & bad_discriminated_taxonomies
    num_overlapping_taxonomies = len(overlapping_taxonomies)

    percentage_good_discrimination = (good_discrimination / total_clusters * 100) if total_clusters > 0 else 0
    num_unique_taxonomies = len(unique_taxonomies)
    num_good_discriminated_taxonomies = len(good_discriminated_taxonomies)

    percentage_good_taxonomies = (num_good_discriminated_taxonomies / num_unique_taxonomies * 100) if num_unique_taxonomies > 0 else 0

    with open(output_stats_file, 'w') as f:
        f.write(f"Number of unique taxonomies: {num_unique_taxonomies}\n")
        f.write(f"Number of well-discriminated taxonomies: {num_good_discriminated_taxonomies}\n")
        f.write(f"Percentage of well-discriminated taxonomies: {percentage_good_taxonomies:.2f}%\n")
        f.write(f"Number of taxonomies both well and poorly discriminated: {num_overlapping_taxonomies}\n")
        f.write(f"Total number of filtered clusters: {total_clusters}\n")
        f.write(f"Number of clusters with good discrimination: {good_discrimination}\n")
        f.write(f"Number of clusters with poor discrimination: {bad_discrimination}\n")
        f.write(f"Percentage of clusters with good discrimination: {percentage_good_discrimination:.2f}%\n\n")

        f.write("Poorly discriminated clusters:\n")
        for cluster, content in bad_clusters:
            f.write(f"Cluster {cluster}:\n")
            for line in content.values():
                f.write(f"{line}\n")
            f.write("\n")

        f.write("Taxonomies both well and poorly discriminated:\n")
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
    print(f"Number of well-discriminated taxonomies: {num_good_discriminated_taxonomies}")
    print(f"Percentage of well-discriminated taxonomies: {percentage_good_taxonomies:.2f}%")
    print(f"Number of taxonomies both well and poorly discriminated: {num_overlapping_taxonomies}")

    return (bad_clusters,
            taxonomy_good_clusters,
            taxonomy_bad_clusters,
            overlapping_taxonomies)

def get_html_stats_content(stats_file):
    """
    Reads a statistics text file and converts it into HTML (without the <html>, <body> tags).
    Returns a string containing HTML for the statistics.
    """
    html_content = []
    with open(stats_file, 'r') as txt_file:
        lines = txt_file.readlines()

        for line in lines:
            line = line.strip()
            # We can make basic paragraphs or headings for each line
            if (line.startswith("Number") or
                line.startswith("Percentage") or
                line.startswith("Total") or
                line.startswith("Poorly discriminated clusters:") or
                line.startswith("Taxonomies both well and poorly discriminated:")):
                html_content.append(f"<p><b>{line}</b></p>\n")
            elif line.startswith("Cluster"):
                html_content.append(f"<h4>{line}</h4>\n")
            elif line:
                # If the line has the format: seq_id \t taxid \t taxonomy
                parts = line.split("\t")
                if len(parts) == 3:
                    seq_id, taxid, taxonomy = parts
                    html_content.append(f"<p>{seq_id} <a href='https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}' target='_blank'>{taxid}</a> {taxonomy}</p>\n")
                else:
                    html_content.append(f"<p>{line}</p>\n")
            else:
                html_content.append("<br>\n")

    return "".join(html_content)

def get_html_cluster_content(cluster_file):
    """
    Reads a cluster file (filtered or rejected) and converts it into HTML list format
    (without the <html>, <body> tags).
    Returns a string containing HTML for the clusters.
    """
    html_content = []
    with open(cluster_file, 'r') as txt_file:
        lines = txt_file.readlines()

        in_cluster = False
        for line in lines:
            line = line.strip()
            if line.startswith("Cluster"):
                # Close previous <ul> if needed
                if in_cluster:
                    html_content.append("</ul>\n")
                html_content.append(f"<h3>{line}</h3>\n<ul>\n")
                in_cluster = True
            else:
                parts = line.split("\t")
                if len(parts) == 3:
                    seq_id, taxid, taxonomy = parts
                    html_content.append(f"<li>{seq_id} <a href='https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}' target='_blank'>{taxid}</a> {taxonomy}</li>\n")

        if in_cluster:
            html_content.append("</ul>\n")

    return "".join(html_content)

def generate_single_html_with_tabs(stats_file, cluster_corrected_file, log_file, output_html_file):
    """
    Generates a single HTML file with 3 tabs:
      1) Statistics
      2) Corrected Clusters
      3) Rejected Clusters
    The content for each tab is built from the provided text files.
    """
    stats_html = get_html_stats_content(stats_file)
    corrected_html = get_html_cluster_content(cluster_corrected_file)

    # If the user provided a log file, we get its content; otherwise, we make an empty string
    if log_file:
        rejected_html = get_html_cluster_content(log_file)
    else:
        rejected_html = "<p>No rejected clusters were logged.</p>"

    # Build the tabbed HTML structure
    html_output = []
    html_output.append("<html>\n<head>\n")
    html_output.append("""
<style>
/* Basic tab styling */
.tab {
  overflow: hidden;
  border-bottom: 1px solid #ccc;
}

.tab button {
  background-color: #f1f1f1;
  border: none;
  outline: none;
  cursor: pointer;
  padding: 10px 20px;
  transition: 0.3s;
  font-size: 16px;
  margin-right: 5px;
}

.tab button:hover {
  background-color: #ddd;
}

.tab button.active {
  background-color: #ccc;
}

.tabcontent {
  display: none;
  padding: 10px;
  border: 1px solid #ccc;
  border-top: none;
}
</style>

<script>
function openTab(evt, tabName) {
  // Hide all tabcontent
  var tabcontents = document.getElementsByClassName("tabcontent");
  for (var i = 0; i < tabcontents.length; i++) {
    tabcontents[i].style.display = "none";
  }
  // Remove the 'active' class from all buttons
  var tablinks = document.getElementsByTagName("button");
  for (var i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }
  // Show the current tab, and add 'active' class to the button that opened the tab
  document.getElementById(tabName).style.display = "block";
  evt.currentTarget.className += " active";
}
</script>
""")

    html_output.append("</head>\n<body>\n")

    # Tab buttons
    html_output.append("""
<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'Statistics')">Statistics</button>
  <button class="tablinks" onclick="openTab(event, 'Corrected')">Corrected Clusters</button>
  <button class="tablinks" onclick="openTab(event, 'Rejected')">Rejected Clusters</button>
</div>
""")

    # Statistics content
    html_output.append(f"""
<div id="Statistics" class="tabcontent">
  <h2>Statistics</h2>
  {stats_html}
</div>
""")

    # Corrected content
    html_output.append(f"""
<div id="Corrected" class="tabcontent">
  <h2>Corrected Clusters</h2>
  {corrected_html}
</div>
""")

    # Rejected content
    html_output.append(f"""
<div id="Rejected" class="tabcontent">
  <h2>Rejected Clusters</h2>
  {rejected_html}
</div>
""")

    # Set the default open tab (Statistics)
    html_output.append("""
<script>
// By default, open the Statistics tab
document.getElementById("Statistics").style.display = "block";
</script>
""")

    html_output.append("</body>\n</html>")

    # Write the final HTML
    with open(output_html_file, 'w') as out_html:
        out_html.write("".join(html_output))

    print(f"Single HTML with tabs generated: {output_html_file}")
    

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################


def main():
    parser = argparse.ArgumentParser(description="Filter clusters with optional cleaning and all_species or selected_species options.")
    parser.add_argument("-c", "--cluster_file", type=str, required=True, help="Input file containing cluster information.")
    parser.add_argument("-o", "--output_stats_file", type=str, required=True, help="Output file for the statistics report.")
    parser.add_argument("-r", "--output_corrected_file", type=str, required=True, help="Output file for the filtered clusters.")
    parser.add_argument("-i", "--include_keywords", type=str, nargs='*', help="Keywords that must be present in taxonomy to keep a line.")
    parser.add_argument("-e", "--exclude_keywords", type=str, nargs='*', help="Keywords that, if present, will cause a line to be excluded.")
    parser.add_argument("-l", "--log_file", type=str, help="Output file to log rejected clusters (optional).")
    parser.add_argument("--clean_words", action='store_true', help="Clean clusters based on a predefined dictionary of words.")
    parser.add_argument("--all_species", action='store_true', help="If set, keeps the entire cluster if it contains at least one valid line.")
    parser.add_argument("--selected_species", action='store_true', help="If set, removes clusters with only excluded lines but keeps clusters with mixed lines, removing only excluded lines.")

    args = parser.parse_args()

    # Force the user to use either --all_species or --selected_species if -i or -e are used
    if (args.include_keywords or args.exclude_keywords) and not (args.all_species or args.selected_species):
        raise ValueError("You must use --all_species or --selected_species when using -i or -e.")

    # User cannot use both all_species and selected_species
    if args.all_species and args.selected_species:
        raise ValueError("You cannot use --all_species and --selected_species at the same time.")

    try:
        filtered_clusters, rejected_clusters, all_taxonomies = filter_lines_by_keywords(
            args.cluster_file,
            include_keywords=args.include_keywords,
            exclude_keywords=args.exclude_keywords,
            clean_words=args.clean_words,
            all_species=args.all_species,
            selected_species=args.selected_species
        )

        # Write the filtered clusters to file
        with open(args.output_corrected_file, 'w') as f:
            for cluster, content in filtered_clusters:
                f.write(f"Cluster {cluster}:\n")
                for line in content.values():
                    f.write(f"{line}\n")

        # Write the rejected clusters to the log file if provided
        if args.log_file:
            with open(args.log_file, 'w') as log:
                for cluster, content in rejected_clusters:
                    log.write(f"Cluster {cluster}:\n")
                    for line in content.values():
                        log.write(f"{line}\n")

        # Calculate statistics and produce the stats file
        bad_clusters, taxonomy_good_clusters, taxonomy_bad_clusters, overlapping_taxonomies = calculate_statistics(
            filtered_clusters,
            args.output_stats_file,
            "uniq_taxo.txt",
            "uniq_taxo_good_discriminated.txt"
        )

        # Now produce a single HTML with 3 tabs:
        #   1. Stats
        #   2. Corrected
        #   3. Rejected
        output_html_file = "final_report.html"
        generate_single_html_with_tabs(
            args.output_stats_file,
            args.output_corrected_file,
            args.log_file,
            output_html_file
        )

        print(f"Statistics written to {args.output_stats_file}")
        print(f"Filtered clusters written to {args.output_corrected_file}")
        if args.log_file:
            print(f"Rejected clusters logged to {args.log_file}")
        print(f"Complete HTML report with tabs: {output_html_file}")

    except ValueError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
