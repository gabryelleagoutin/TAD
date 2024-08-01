import os
import argparse
from Bio.Seq import Seq
from Bio import SeqIO

def parse_primer_info(line):
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
    return filename.split('_')[0]

def get_alignment_size(og_id, alignment_folder):
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

def has_single_gc_clamp(sequence):
    '''Function to check if the sequence has a single GC clamp'''
    last_six_bases = sequence[-6:]
    return (last_six_bases.count('G') + last_six_bases.count('C')) == 1

def amplicon_score(amplicon_size, amplicon_min_size):
    return round(((amplicon_size - amplicon_min_size) / 22), 2)
    
def process_files(input_files, alignment_folder, amplicon_min_size, amplicon_max_size):
    all_primer_pairs = []
    
    for tsv_file in input_files:
        output_file = tsv_file.replace('.tsv', '_couple.tsv')  # Génère le nom du fichier de sortie
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
                gc_clamp_RC_B = has_single_gc_clamp(reverse_complement_B)
                gc_last_thirty_percent_RC_B = count_gc_in_last_thirty_percent(reverse_complement_B)
                primer_pair = '\t'.join(
                    [primers[i]['OG_ID'], primers[i]['NumberOfSeq'], primers[i]['SpeciesCount'], primers[i]['PercentSingleCopy'],
                     primers[i]['GeneName'], str(alignment_size)] + primer1_info.split('\t')[5:] +
                    primer2_info[5:] + [reverse_complement_B, str(gc_clamp_RC_B), str(gc_last_thirty_percent_RC_B), str(potential_amplicon_size), str(amplicon_size_score), str(total_score)])
                all_primer_pairs.append(primer_pair + '\n')

    if all_primer_pairs:
        header = "OG_ID\tNumberOfSeq\tSpeciesCount\tPercentSingleCopy\tGeneName\tAlignement_size\tPrimer_A\tPosition_A\tPrimer_Size_A\tNumber_matching_A\tPercentage_NM_A\tScore_Percentage_NM_A\tDegenerescence_A\tTm_A_max\tTm_A_min\tGC_percentage_fraction_A\tGC_percentage_max_A\tGC_percentage_min_A\tGC_in_last_thirty_percent_A\tEnds_with_T_A\tSelf_Complementarity_A\tGC_clamp_A\tPrimer_B\tPosition_B\tPrimer_Size_B\tNumber_matching_B\tPercentage_NM_B\tScore_Percentage_NM_B\tDegenerescence_B\tTm_B_max\tTm_B_min\tGC_percentage_fraction_B\tGC_percentage_max_B\tGC_percentage_min_B\tGC_in_last_thirty_percent_B\tEnds_with_T_B\tSelf_Complementarity_B\tGC_clamp2\tReverse_Complement_B\tGC_clamp_RC_B\tGC_last_6_bases_RC_B\tpotential_amplicon_size\tAmplicon_score\tTotal_score"
        with open(output_file, 'w') as out_file:
            out_file.write(header + '\n')
            for pair in all_primer_pairs:
                out_file.write(pair)
    else:
        open(output_file, 'w').close()  # Create an empty file if no primer pairs found

def main():
    parser = argparse.ArgumentParser(description='Process some TSV files.')
    parser.add_argument("-i", "--input_files", nargs='+', help="Paths to the input TSV files")
    parser.add_argument('-f', '--alignment_folder', type=str, required=True, help='The folder containing alignment files')
    parser.add_argument('--amplicon_min_size', type=int, default=50, help='Minimum size of the amplicon')
    parser.add_argument('--amplicon_max_size', type=int, default=600, help='Maximum size of the amplicon')
    args = parser.parse_args()

    process_files(args.input_files, args.alignment_folder, args.amplicon_min_size, args.amplicon_max_size)

if __name__ == "__main__":
    main()
