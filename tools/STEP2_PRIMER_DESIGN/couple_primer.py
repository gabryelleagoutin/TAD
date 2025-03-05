#!/usr/bin/env python

import os
import argparse
from Bio.Seq import Seq
from Bio import SeqIO

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

def parse_primer_info(line):
    '''Parse a line from the primer information file and return a dictionary with primer details.'''
    primer_info = line.split('\t')
    primer = {
        'OG_ID': primer_info[0],
        'NumberOfSeq': primer_info[1],
        'SpeciesCount': primer_info[2],
        'PercentSingleCopy': primer_info[3],
        'GeneName': primer_info[4],
        'Primer': primer_info[5],
        'Position': int(primer_info[6]),
        'Score_Percentage_NM': primer_info[10],
        'Info': line.strip()
    }
    return primer

def get_og_id(filename):
    '''Extract the OG_ID from the filename by splitting at the underscore.'''
    return filename.split('_')[0]

def get_alignment_size(og_id, alignment_folder):
    '''Get the size of sequences in the alignment files for a given OG_ID.'''
    for file_name in os.listdir(alignment_folder):
        if file_name.startswith(og_id):
            sequence_path = os.path.join(alignment_folder, file_name)
            with open(sequence_path) as fasta_file:
                records = SeqIO.parse(fasta_file, 'fasta')
                sizes = set()
                for record in records:
                    sizes.add(len(record.seq))
                if len(sizes) == 1:
                    return sizes.pop()
                else:
                    return "Different sizes for sequences"
    return None

def count_gc_in_last_thirty_percent(sequence):
    '''Function to count GC bases in the last thirty percent of the sequence'''
    thirty_percent_length = int(len(sequence) * 0.3)
    last_thirty_percent_bases = sequence[-thirty_percent_length:]
    gc_count = last_thirty_percent_bases.count('G') + last_thirty_percent_bases.count('C')
    
    return gc_count

def amplicon_score(amplicon_size, amplicon_min_size):
    '''Calculate a score for the amplicon size based on the size.'''
    return round(((amplicon_size - amplicon_min_size) / 22), 2)
    
def process_files(input_files, alignment_folder, amplicon_min_size, amplicon_max_size):
    '''Process each input TSV file, find primer pairs, and save the results to output files.'''
    all_primer_pairs = []
    
    for tsv_file in input_files:
        output_file = tsv_file.replace('.tsv', '_couple.tsv') 
        with open(tsv_file, 'r') as file:
            lines = file.readlines()

        primers = [parse_primer_info(line) for line in lines[1:]]  # Skip the header line
        primers.sort(key=lambda x: x['Position'])

        for i in range(len(primers)):
            for j in range(i + 1, len(primers)):
                primer1_info = primers[i]['Info']
                primer2_info = primers[j]['Info'].split('\t')
                position_1 = int(primer1_info.split('\t')[6])
                position_2 = int(primer2_info[6])
                primer_length_1 = len(primer1_info.split('\t')[5])
                potential_amplicon_size = position_2 - position_1 - primer_length_1
                
                if potential_amplicon_size < amplicon_min_size or potential_amplicon_size > amplicon_max_size:
                    continue  # Skip this pair of primers if it doesn't meet the size criteria
                    
                amplicon_size_score = amplicon_score(potential_amplicon_size, amplicon_min_size)

                score_percentage_nm_a = float(primers[i]['Score_Percentage_NM'])
                score_percentage_nm_b = float(primer2_info[10])
                min_score_percentage_nm = min(score_percentage_nm_a, score_percentage_nm_b)
                total_score = round(min_score_percentage_nm + amplicon_size_score, 2)

                og_id = get_og_id(primers[i]['OG_ID'])
                alignment_size = get_alignment_size(og_id, alignment_folder)
                reverse_complement_B = str(Seq(primer2_info[5]).reverse_complement())
                gc_last_thirty_percent_RC_B = count_gc_in_last_thirty_percent(reverse_complement_B)
                primer_pair = '\t'.join(
                    [primers[i]['OG_ID'], primers[i]['NumberOfSeq'], primers[i]['SpeciesCount'], primers[i]['PercentSingleCopy'],
                     primers[i]['GeneName'], str(alignment_size)] + primer1_info.split('\t')[5:] +
                    primer2_info[5:] + [reverse_complement_B, str(gc_last_thirty_percent_RC_B), str(potential_amplicon_size), str(amplicon_size_score), str(total_score)])
                all_primer_pairs.append(primer_pair + '\n')

    if all_primer_pairs:
        header = "OG_ID\tNumberOfSeq\tSpeciesCount\tPercentSingleCopy\tGeneName\tAlignement_size\tPrimer_A\tPosition_A\tPrimer_Size_A\tNumber_matching_A\tPercentage_NM_A\tScore_Percentage_NM_A\tDegenerescence_A\tTm_A_max\tTm_A_min\tGC_percentage_fraction_A\tGC_percentage_max_A\tGC_percentage_min_A\tGC_in_last_thirty_percent_A\tEnds_with_T_A\tSelf_Complementarity_A\tGC_clamp_A\tPrimer_B\tPosition_B\tPrimer_Size_B\tNumber_matching_B\tPercentage_NM_B\tScore_Percentage_NM_B\tDegenerescence_B\tTm_B_max\tTm_B_min\tGC_percentage_fraction_B\tGC_percentage_max_B\tGC_percentage_min_B\tGC_in_last_thirty_percent_B\tEnds_with_T_B\tSelf_Complementarity_B\tGC_clamp2\tReverse_Complement_B\tGC_last_trhity_percent_RC_B\tpotential_amplicon_size\tAmplicon_score\tTotal_score"
        with open(output_file, 'w') as out_file:
            out_file.write(header + '\n')
            for pair in all_primer_pairs:
                out_file.write(pair)
    else:
        open(output_file, 'w').close()  # Create an empty file if no primer pairs found
##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

def main():
    parser = argparse.ArgumentParser("""description='realise all possible pairs of primers.
Columns in the output file:
    - $1 OG_ID: ID of the Orthologous Group
    - $2 NumberOfSeq: Total number of sequences in the fasta file
    - $3 SpeciesCount: Number of different species in this OG
    - $4 PercentSingleCopy: Percentage of proteins that are single-copy in the OG
    - $5 GeneName
    - $6 Alignement_size : size of the aligned sequences, gives an estimate of the size of the sequences.
    - $7 Primer_A : DNA sequence of the forward primer    
    - $8 Position_A: Forward primer position on sequence alignment
    - $9 Primer_Size_A: Forward primer size
    - $10 Number_matching_A : Number of sequences matching this primer
    - $11 Percentage_NM : Percentage. Number of sequences caught by the primer as a function of the total number of sequences
    - $12 Score_Percentage_NM: Score calculated as the difference between the percentage observed and a predefined threshold.
    - $13 Degenerescence_A: Degeneracy of the primer      
    - $14 Tm_A_max: 2*(A+T)+4*(G+C) with max of GC        
    - $15 Tm_A_min: 2*(A+T)+4*(G+C) with max of AT     
    - $16 GC_percentage_fraction_A: Fractional GC content of the primer with degeneracy
    - $17 GC_percentage_max_A: Maximum GC content of the primer     
    - $18 GC_percentage_min_A: Minimum GC content of the primer     
    - $19 GC_in_last_thirty_percent_A: Number of GC bases in the last thirty percent base of the primer.To check GC content at the end of the primer     
    - $20 Ends_with_T_A: Whether the primer ends with 'T'   
    - $21 Self_Complementarity_A: Whether the primer is self-complementary  
    - $22 GC_clamp_A : Whether the primer has a single GC clamp     
    - $23 Primer_B: DNA sequence of the forward primer       
    - $24 Position_B      
    - $25 Primer_Size_B
    - $26 Number_matching_B
    - $27 Percentage_NM_B
    - $28 Score_Percentage_NM_B
    - $29 Degenerescence_B        
    - $30 Tm_B_max        
    - $31 Tm_B_min        
    - $32 GC_percentage_fraction_B        
    - $33 GC_percentage_max_B     
    - $34 GC_percentage_min_B     
    - $35 GC_in_last_thirty_percent_B     
    - $36 Ends_with_T_B
    - $37 Self_Complementarity_B   
    - $38 GC_clamp2       
    - $39 Reverse_Complement_B: reverse complement primer for in silico pcr 
    - $40 GC_last_trhity_percent_RC_B 
    - $41 potential_amplicon_size: amplicon size between the two primers
    - $42 Amplicon_score: +1 every 22 bases (amplicon_size - amplicon_min_size) / 22
    - $43 Total_score: Score_Percentage_NM the lowest score of Score_Percentage_NM + the amplicon_score. We take the weakest base, because that's the one that would catch the most primers.""", formatter_class=argparse.RawTextHelpFormatter,
    epilog="python couple_primer.py -i $file -f alignment/ --amplicon_min_size 150 --amplicon_max_size 590")
    parser.add_argument("-i", "--input_files", nargs='+', help="Paths to the input TSV files")
    parser.add_argument('-f', '--alignment_folder', type=str, required=True, help='The folder containing alignment files')
    parser.add_argument('--amplicon_min_size', type=int, default=150, help='Minimum size of the amplicon')
    parser.add_argument('--amplicon_max_size', type=int, default=490, help='Maximum size of the amplicon')
    args = parser.parse_args()

    process_files(args.input_files, args.alignment_folder, args.amplicon_min_size, args.amplicon_max_size)

if __name__ == "__main__":
    main()
