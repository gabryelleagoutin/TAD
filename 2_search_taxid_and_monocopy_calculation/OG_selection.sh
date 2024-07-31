#!/bin/bash

# Function to display help
display_help() {
    echo "Usage: $0 -p <percent_single_copy> [-c <TargetSpecies_Count>] [-t <TargetSpecies_Percentage>] -o <output_file> <data_file>"
    echo "Options:"
    echo "  -p    specifies the threshold X for percent_single_copy (mandatory): genes are in single copy in at least X% of species"
    echo "  -c    specifies the threshold for TargetSpecies_Count (optional): If you know the number of species present in the OrthoDB data, you can provide a minimum. You can likely find this number in the odb11v0_levels.tab file by grepping for the name of your taxonomic rank. The last column gives you the number of species associated with it."
    echo "  -t    specifies the threshold for TargetSpecies_Percentage (optional): In this ortholog group, there are X% of my species of interest. The rest are not of interest."
    echo "  -o    specifies the output file name (mandatory)"
    echo "  -h    displays this help message"
    exit 0
}

# Initialize default values
percent_single_copy=""
TargetSpecies_Count=""
TargetSpecies_Percentage=""
output_file=""

# Check command line options
while getopts ":p:c:t:o:h" opt; do
    case $opt in
        p)
            percent_single_copy="$OPTARG"
            ;;
        c)
            TargetSpecies_Count="$OPTARG"
            ;;
        t)
            TargetSpecies_Percentage="$OPTARG"
            ;;
        o)
            output_file="$OPTARG"
            ;;
        h)
            display_help
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

# Check if percent_single_copy and output_file have been specified (mandatory)
if [ -z "$percent_single_copy" ] || [ -z "$output_file" ]; then
    echo "The arguments -p (percent_single_copy) and -o (output_file) are mandatory."
    exit 1
fi

# Ignore options processed by getopts
shift $((OPTIND-1))

# Check if the number of positional arguments is correct
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 -p <percent_single_copy> [-c <TargetSpecies_Count>] [-t <TargetSpecies_Percentage>] -o <output_file> <data_file>"
    exit 1
fi

# Store the file name as a variable
data_file="$1"

# Check if the file exists
if [ ! -e "$data_file" ]; then
    echo "The file $data_file does not exist."
    exit 1
fi

# Use Awk to filter lines based on the given criteria
awk -F '\t' -v col5="$percent_single_copy" -v col10="$TargetSpecies_Count" -v col11="$TargetSpecies_Percentage" '$5 >= col5 && $10 >= col10 && $11 >= col11' "$data_file" > "$output_file"
