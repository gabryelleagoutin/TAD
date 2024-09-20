#!/usr/bin/env python
import argparse
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt, GC123
import csv
import os

__author__ = 'Gabryelle Agoutin - INRAE'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'gabryelle.agoutin@inrae.fr'
__status__ = 'prod'

def validate_file_exists(file_path):
    '''Function to validate if a file exists'''
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Le fichier {file_path} n'existe pas.")

def process_file_tsv(file_path):
    '''Function to process a TSV file and yield its content'''
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        header = next(reader)  # Ignorer l'en-tête
        for line in reader:
            yield line

def calculate_gc_percentage_with_degeneracy(sequence):
    '''Function to calculate GC percentage considering base degeneracy'''
    base_deg = {'G': 1, 'C': 1, 'S': 1, 'R': 0.5, 'Y': 0.5, 'K': 0.5, 'M': 0.5, 'B': 0.667, 'V': 0.667, 'D': 0.333, 'H': 0.333, 'N': 0.5}

    gc_count = sum(base_deg.get(base, 0) for base in sequence)
    total_bases = len(sequence)

    gc_percentage = gc_count / total_bases if total_bases > 0 else 0.0
    return round((gc_percentage)*100, 2)

def calculate_gc_percentage_max(sequence):
    '''Function to calculate maximum GC percentage'''
    allowed_bases = {'G', 'C', 'R', 'Y', 'S', 'K', 'M', 'B', 'V', 'D', 'H', 'N'}

    gc_count = sum(1 for base in sequence if base in allowed_bases)
    total_bases = len(sequence)

    gc_percentage = gc_count / total_bases if total_bases > 0 else 0.0
    return round((gc_percentage)*100, 2)

def calculate_gc_percentage_min(sequence):
    '''Function to calculate minimum GC percentage'''
    gc_count = sum(1 for base in sequence if base in {'G', 'C', 'S'})
    total_bases = len(sequence)

    gc_percentage = gc_count / total_bases if total_bases > 0 else 0.0
    return round((gc_percentage)*100, 2)

def percent_NM(number_matching, number_of_seq):
    '''percentage matching sequence'''
    
    if number_of_seq > 0:
        nm_percentage = (number_matching / number_of_seq) * 100
        return round(nm_percentage, 2)
    else:
        return 0.0

def score_percentage_nm(nm_percentage, threshold):
    if nm_percentage < threshold:
        return "under_threshold"
    else:
        return round((nm_percentage - threshold), 2)

def count_gc_in_last_thirty_percent(sequence):
    '''Function to count GC bases in the last thirty percent of the sequence'''
    thirty_percent_length = int(len(sequence) * 0.3)
    last_thirty_percent_bases = sequence[-thirty_percent_length:]
    gc_count = last_thirty_percent_bases.count('G') + last_thirty_percent_bases.count('C')

    return gc_count

def has_single_gc_clamp(sequence, num_bases=5):
    """
    Function to check if there is at least one G, C, or S in the last 'num_bases' bases.
    
    Parameters:
    - sequence (str): DNA sequence that may contain degenerate bases.
    - num_bases (int): The number of bases to check from the end of the sequence (default is 5).
    
    Returns:
    - bool: True if there is a GC clamp (at least one G, C, or S), False otherwise.
    """
    last_bases = sequence[-num_bases:]
    return any(base in ['G', 'C', 'S'] for base in last_bases)
    
    
def ends_with_t(sequence):
    '''Function to check if the sequence ends with T'''
    return sequence[-1] == 'T'

def self_complementarity(sequence):
    '''Function to check if the sequence is self-complementary'''
    return sequence == sequence.reverse_complement()

def Tm_maxi(sequence):
    tm = 0 
    nuc_GC = 0
    nuc_AT = 0
    bases_GC = {'G', 'C', 'R', 'Y', 'S', 'K', 'M', 'B', 'V', 'D', 'H', 'N'}
    for base in sequence:
        if base in bases_GC:
            nuc_GC += 1
        else:
            nuc_AT +=1        
    tm = float(2*(nuc_AT)+4*(nuc_GC))
    return(tm)
    
def Tm_mini(sequence):
    tm = 0 
    nuc_GC = 0
    nuc_AT = 0
    bases_AT = {'A', 'T','W' 'R', 'Y', 'K', 'M', 'B', 'V', 'D', 'H', 'N'}
    for base in sequence:
        if base in bases_AT:
            nuc_AT += 1
        else:
            nuc_GC +=1        
    tm = float(2*(nuc_AT)+4*(nuc_GC))
    return(tm)

def process_line(columns, og_id, og_info, nm_threshold, tm_max_threshold, tm_min_threshold):
    '''Function to process a line and extract features'''
    position = columns[0]
    primer = Seq(columns[5])
    primer_size = len(primer)
    number_matching = int(columns[6])
    degenerescence = int(columns[4])
    tm_max = Tm_maxi(primer)
    tm_min = Tm_mini(primer)
    gc_percentage = calculate_gc_percentage_with_degeneracy(primer)
    gc_percentage_max = calculate_gc_percentage_max(primer)
    gc_percentage_min = calculate_gc_percentage_min(primer)
    gc_in_last_thirty_percent = count_gc_in_last_thirty_percent(primer)
    ends_with_t_flag = ends_with_t(primer)
    self_complementarity_flag = self_complementarity(primer)
    gc_clamp_flag = has_single_gc_clamp(primer)
    og_id= str(og_id) # sinon j'ai erreur : invalid literal for int() with base 10: '72971at1578'
    og_data = og_info.get(og_id, {})
    number_of_seq = og_data.get("NumberOfSeq", "")
    percentage_nm = percent_NM(int(number_matching), int(number_of_seq))
    score_percentage = score_percentage_nm(float(percentage_nm), nm_threshold)
    
    if score_percentage == "under_threshold" or tm_max > tm_max_threshold or tm_min < tm_min_threshold:
        return None

    return {
        "OG_ID": og_id,
        "NumberOfSeq": og_data.get("NumberOfSeq", ""),
        "SpeciesCount": og_data.get("SpeciesCount", ""),
        "PercentSingleCopy": og_data.get("percent_single_copy", ""),
        "GeneName": og_data.get("gene_name", ""),
        "Primer": str(primer),
        "Position": position,
        "Primer_Size": primer_size,
        "Number_matching": number_matching,
        "Percentage_NM": percentage_nm,
        "Score_Percentage_NM": score_percentage,
        "Degenerescence": degenerescence,
        "Tm_max": tm_max,
        "Tm_min": tm_min,
        "GC_percentage_fraction": gc_percentage,
        "GC_percentage_max": gc_percentage_max,
        "GC_percentage_min": gc_percentage_min,
        "GC_in_last_thirty_percent": gc_in_last_thirty_percent,
        "Ends_with_T": ends_with_t_flag,
        "Self_Complementarity": self_complementarity_flag,
        "GC_clamp": gc_clamp_flag
    }

def process_file(file_path, og_id, og_info, nm_threshold, tm_max_threshold, tm_min_threshold):
    '''Function to process a TSV file and return a list of processed lines'''
    results = []
    for line in process_file_tsv(file_path):
        result = process_line(line, og_id, og_info, nm_threshold, tm_max_threshold, tm_min_threshold)
        if result:
            results.append(result)
    return results

def write_output_table(results, output_file):
    '''Function to write the output table to a file in the current directory'''
    current_directory = os.getcwd()
    output_path = os.path.join(current_directory, output_file)
    with open(output_path, 'w') as out_file:
        out_file.write("OG_ID\tNumberOfSeq\tSpeciesCount\tPercentSingleCopy\tGeneName\tPrimer\tPosition\tPrimer_Size\tNumber_matching\tPercentage_NM\tScore_Percentage_NM\tDegenerescence\tTm_max\tTm_min\tGC_percentage_fraction\tGC_percentage_max\tGC_percentage_min\tGC_in_last_thirty_percent\tEnds_with_T\tSelf_Complementarity\tGC_clamp\n")
        for result in results:
            out_file.write("\t".join(str(result[key]) for key in ["OG_ID","NumberOfSeq", "SpeciesCount", "PercentSingleCopy", "GeneName", "Primer", "Position", "Primer_Size", "Number_matching","Percentage_NM", "Score_Percentage_NM", "Degenerescence", "Tm_max", "Tm_min", "GC_percentage_fraction","GC_percentage_max", "GC_percentage_min","GC_in_last_thirty_percent", "Ends_with_T", "Self_Complementarity", "GC_clamp"]) + "\n")

def read_og_info(og_file):
    '''Function to read OG info from a file and store it in a dictionary'''
    og_info = {}
    with open(og_file, 'r') as og_file:
        next(og_file)  # Skip header
        for line in og_file:
            columns = line.strip().split('\t')
            og_id = columns[0]
            species_count = columns[9]
            percent_single_copy = columns[4]
            gene_name = columns[8]
            number_of_seq = columns[11]
            og_info[og_id] = {
                "NumberOfSeq": number_of_seq,
                "SpeciesCount": species_count,
                "percent_single_copy": percent_single_copy,
                "gene_name": gene_name
            }
    return og_info

def main():
    parser = argparse.ArgumentParser(description="""This script processes the degeprime results, extracts various characteristics, and calculates others.

Columns in the output file:
    - $1 OG_ID: ID of the Orthologous Group
    - $2 NumberOfSeq: Total number of sequences in the fasta file
    - $3 SpeciesCount: Number of different species in this OG
    - $4 PercentSingleCopy : Percentage of proteins that are single-copy in the OG
    - $5 Gene_name
    - $6 Primer: DNA sequence of the primer
    - $7 Position : primer position on sequence alignment
    - $8 Primer_Size: Size of the primer
    - $9 Number_matching: Number of sequences matching this primer
    - $10 Percentage_NM : Percentage. Number of sequences caught by the primer as a function of the total number of sequences
    - $11 Score_Percentage_NM : Score calculated as the difference between the percentage observed and a predefined threshold. 
    - $12 Degenerescence: Degeneracy of the primer
    - $13 Tm_max: 2*(A+T)+4*(G+C) with max of GC
    - $14 Tm_min: 2*(A+T)+4*(G+C) with max of AT
    - $15 GC_percentage_fraction: Fractional GC content of the primer with degeneracy
    - $16 GC_percentage_max: Maximum GC content of the primer
    - $17 GC_percentage_min: Minimum GC content of the primer
    - $18 GC_in_last_thirty_percent: Number of GC bases in the last thirty percent base of the primer.To check GC content at the end of the primer
    - $19 Ends_with_T: Whether the primer ends with 'T'
    - $20 Self_Complementarity: Whether the primer is self-complementary
    - $21 GC_clamp: Whether the primer has a single GC clamp""", formatter_class=argparse.RawTextHelpFormatter,
    epilog="Exemple: python process_primers_stat.py -i degeprime_result/concatenated_* -og ../2_search_taxid_and_monocopy_calculation/OG_selected_1578.tab -o result_stat_primers -nm 80 -tm_max 70 -tm_min 50")
    parser.add_argument("-i", "--input_files", nargs='+', help="Paths to the input TSV files")
    parser.add_argument("-og", "--og_file", help="Path to the OG information file")
    parser.add_argument("-o", "--output_dir", help="Répertoire de sortie pour les fichiers de sortie")
    parser.add_argument("-nm", "--nm_threshold", type=float, default=80, help="Threshold for percentage NM filtering")
    parser.add_argument("-tm_max", "--tm_max_threshold", type=float, default=70, help="Threshold for maximum temperature filtering")
    parser.add_argument("-tm_min", "--tm_min_threshold", type=float, default=50, help="Threshold for minimum temperature filtering")
    args = parser.parse_args()

    og_info = read_og_info(args.og_file)

    for input_file in args.input_files:
        validate_file_exists(input_file)

        try:
            # Extract OG ID from filename
            og_id = os.path.basename(input_file).split('_')[1].split('.')[0]

            results = process_file(input_file, og_id, og_info, args.nm_threshold, args.tm_max_threshold, args.tm_min_threshold)
            
            # Create output file name based on input file name
            output_file = os.path.join(args.output_dir, os.path.basename(os.path.splitext(input_file)[0]) + "_stat_primer.tsv")
            # Write the results to the output file  
            write_output_table(results, output_file)
            print(f"Results written to {output_file}")
        except Exception as e:
            print(f"An error occurred while processing {input_file}: {e}")

if __name__ == "__main__":
    main()
