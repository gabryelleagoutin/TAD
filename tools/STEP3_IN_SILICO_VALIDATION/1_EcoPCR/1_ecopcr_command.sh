#!/bin/bash

# Default variables
s_array_file='ecopcr_commands.sarray'
ecopcr_outdir='result'

# Function to display script usage
usage() {
    echo "Usage: $0 -a '<primer list>' -d '<path to .adx files>'"
    echo "  -a   Comma-separated list of primer pairs, each consisting of a forward primer and a reverse primer. Primers are separated by a space and surrounded by an inverted comma. e.g. ‘TSRTCAAGAACRTBGARR SAYGTYCTGBACRTCRTC,TSRTCAAGAACRTBGARR AYGTYCTGNACRTCGTCV’."
    echo "  -d   Access path to the database for EcoPCR, the folder containing the .adx files. e.g. /BD_TaxonMarker/BD_ecoPCR/ecoPCR_db/assemblies/"
    exit 1
}

# Parse arguments
while getopts "a:d:" opt; do
    case ${opt} in
        a)
            primer_list=$OPTARG
            ;;
        d)
            ecopcr_assembly_dir=$OPTARG
            ;;
        *)
            usage
            ;;
    esac
done

# Ensure both arguments are provided
if [ -z "$primer_list" ] || [ -z "$ecopcr_assembly_dir" ]; then
    usage
fi

mkdir -p $ecopcr_outdir

# Split the primer list by commas
IFS=',' read -r -a primers <<< "$primer_list"

# Initialize a counter for numbering directories
counter=1

# Loop through each primer pair
for primer_pair in "${primers[@]}"
do
    # Split into forward and reverse primers
    IFS=' ' read -r primer_f primer_r <<< "$primer_pair"

    # Name the output directory with a number and the primers
    outdir="${counter}_${primer_f}-${primer_r}"
    mkdir -p "$ecopcr_outdir/$outdir"

    # Loop through each assembly file (.adx)
    for assembly in $ecopcr_assembly_dir/*.adx
    do
        name=$(basename "$assembly" | sed 's/\.adx$//')

        # Add the ecoPCR command to the s_array_file
        echo "module load bioinfo/ecoPCR/1.0.1; ecoPCR -d $(echo "$assembly" | sed 's/\.adx$//') $primer_f $primer_r > $ecopcr_outdir/$outdir/${name}.ecopcr;" >> $s_array_file
    done

    # Increment the counter after processing the primer pair for all assemblies
    ((counter++))

done
