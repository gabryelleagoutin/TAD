#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

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


def load_data(file_name):
    """Load data from a text file into a set.

    Args:
        file_name (str): The path to the text file.

    Returns:
        set: A set containing the unique lines from the file.
    """
    with open(file_name, 'r') as file:
        return set(line.strip() for line in file)

def save_set_to_file(output_file, data_set):
    """Save a set of data to a text file with a title.

    Args:
        output_file (str): The path to the output text file.
        title (str): Title to be added at the beginning of the file.
        data_set (set): The set of data to be saved.
    """
    with open(output_file, 'w') as file:
        for item in data_set:
            file.write(f"{item}\n")
    print(f"Data saved to '{output_file}'")
    
    
def create_venn_diagrams(data, output_file):
    """Create and save Venn diagrams with descriptions.

    Args:
        data (dict): A dictionary with file names as keys and sets of data as values.
        output_file (str): The path to the output image file.
    """
    # Create figure and subplots
    fig, axs = plt.subplots(1, 2, figsize=(16, 8))  # 1 row, 2 columns
    fig.subplots_adjust(wspace=0.4, left=0.05, right=0.95, top=0.85, bottom=0.2)  # Adjust margins and spacing

    # Diagram for uniq_taxo.txt and uniq_taxo_amplicon.txt
    venn2([data['uniq_taxo.txt'], data['uniq_taxo_amplicon.txt']],
          set_labels=('target gene', 'amplicon gene'),
          ax=axs[0])
    axs[0].set_title('Caught species comparison')
    
    # Description below the first diagram
    axs[0].text(0.5, -0.2, 
                'This diagram compares the species caught by the target gene with those detected by the amplicon gene.\n'
                'Red represents species detected only by the target gene.\n'
                'Green represents species detected only by the amplicon gene.\n'
                'Yellow represents species detected by both genes.',
                ha='center', va='top', fontsize=10, transform=axs[0].transAxes, bbox=dict(facecolor='white', alpha=0.7))

    # Diagram for uniq_taxo_good_discriminated.txt and uniq_taxo_good_discriminated_amplicon.txt
    venn2([data['uniq_taxo_good_discriminated.txt'], data['uniq_taxo_good_discriminated_amplicon.txt']],
          set_labels=('target gene', 'amplicon gene'),
          ax=axs[1])
    axs[1].set_title('Species correctly discriminated')

    # Description below the second diagram
    axs[1].text(0.5, -0.2, 
                'This diagram shows the species that were correctly discriminated by the target gene and the amplicon gene.\n'
                'Red represents species correctly discriminated only by the target gene.\n'
                'Green represents species correctly discriminated only by the amplicon gene.\n'
                'Yellow represents species correctly discriminated by both genes.',
                ha='center', va='top', fontsize=10, transform=axs[1].transAxes, bbox=dict(facecolor='white', alpha=0.7))

    # Save the figure to a file
    plt.savefig(output_file, bbox_inches='tight')  # Use bbox_inches='tight' to fit the layout
    print(f"Venn diagrams with descriptions saved to '{output_file}'")

def save_venn_sets(data):
    """Save the unique and intersecting elements of the Venn diagram to text files.

    Args:
        data (dict): A dictionary with file names as keys and sets of data as values.
        output_prefix (str): Prefix for the output files.
    """
    # Caught species comparison
    only_target_caught = data['uniq_taxo.txt'] - data['uniq_taxo_amplicon.txt']
    only_amplicon_caught = data['uniq_taxo_amplicon.txt'] - data['uniq_taxo.txt']
    intersection_caught = data['uniq_taxo.txt'] & data['uniq_taxo_amplicon.txt']
    
    # Save caught species comparison results
    save_set_to_file('list_only_target_gene_caught.txt', only_target_caught)
    save_set_to_file('list_only_amplicon_gene_caught.txt', only_amplicon_caught)
    save_set_to_file('list_intersection_caught.txt', intersection_caught)

    # Species correctly discriminated
    only_target_discriminated = data['uniq_taxo_good_discriminated.txt'] - data['uniq_taxo_good_discriminated_amplicon.txt']
    only_amplicon_discriminated = data['uniq_taxo_good_discriminated_amplicon.txt'] - data['uniq_taxo_good_discriminated.txt']
    intersection_discriminated = data['uniq_taxo_good_discriminated.txt'] & data['uniq_taxo_good_discriminated_amplicon.txt']

    # Save species correctly discriminated results
    save_set_to_file('list_only_target_gene_discriminated.txt', only_target_discriminated)
    save_set_to_file('list_only_amplicon_gene_discriminated.txt', only_amplicon_discriminated)
    save_set_to_file('list_intersection_discriminated.txt', intersection_discriminated)

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

def main():
    parser = argparse.ArgumentParser(
        description='Create Venn diagrams comparing species data.',
        epilog='Example: python script_name.py --taxo_target_gene path/to/uniq_taxo.txt --taxo_amplicon_gene path/to/uniq_taxo_amplicon.txt --good_discrimination_target_gene path/to/uniq_taxo_good_discriminated.txt --good_discrimination_amplicon path/to/uniq_taxo_good_discriminated_amplicon.txt --output path/to/output_image.png'
    )
    
    # Define arguments
    parser.add_argument('--taxo_target_gene', required=True, help='Path to the uniq_taxo.txt file')
    parser.add_argument('--taxo_amplicon_gene', required=True, help='Path to the uniq_taxo_amplicon.txt file')
    parser.add_argument('--good_discrimination_target_gene', required=True, help='Path to the uniq_taxo_good_discriminated.txt file')
    parser.add_argument('--good_discrimination_amplicon', required=True, help='Path to the uniq_taxo_good_discriminated_amplicon.txt file')
    parser.add_argument('--output', required=True, help='Path to the output image file')

    # Parse arguments
    args = parser.parse_args()

    # Load data from files
    data = {
        'uniq_taxo.txt': load_data(args.taxo_target_gene),
        'uniq_taxo_amplicon.txt': load_data(args.taxo_amplicon_gene),
        'uniq_taxo_good_discriminated.txt': load_data(args.good_discrimination_target_gene),
        'uniq_taxo_good_discriminated_amplicon.txt': load_data(args.good_discrimination_amplicon)
    }

    # Create Venn diagrams
    create_venn_diagrams(data, args.output)
    
    # Save unique and intersecting sets to text files
    save_venn_sets(data)

if __name__ == '__main__':
    main()
