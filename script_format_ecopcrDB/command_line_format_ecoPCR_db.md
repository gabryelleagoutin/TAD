# Retrieve the assembly summary
cp  /work/bank2/users_banks/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/assembly_summary.txt .

# Complete with the header 
cat header.txt
assembly_accession      bioproject      biosample       wgs_master      refseq_category taxid   species_taxid   organism_name   infraspecific_name      isolate version_status  assembly_level  release_type genome_rep       seq_rel_date    asm_name        asm_submitter   gbrs_paired_asm paired_asm_comp ftp_path        excluded_from_refseq    relation_to_type_material       asm_not_live_date       assembly_typegroup    genome_size     genome_size_ungapped    gc_percent      replicon_count  scaffold_count  annotation_provider     annotation_name annotation_date total_gene_count        protein_coding_gene_count    non_coding_gene_count    pubmed_id

cat header.txt assembly_summary.txt> assembly_summary_with_header.txt

# separation into chunks to go faster because of the huge genome
nb_chunk=500
mkdir assembly_chunks
cut -f1 assembly_summary_with_header.txt > ref_and_rep_assemblies_acc_list.txt
split -a 3 -d --number=l/$nb_chunk ref_and_rep_assemblies_acc_list.txt assembly_chunks/ref_and_rep_assemblies.acc_list_chunk.

# the sarray command to format the database
out_db='ecoPCR_db/assemblies'
mkdir -p $out_db
echo '# build ecopcr db on assembly chunk' > build_assemblies_ecopcr_db.sarray  
for chunk_file in assembly_chunks/*chunk*;
do
    command=''
    # parsing chunk file:
    # chunk file has to follow this pattern: <assembly_selection_name>.acc_list_chunk.<chunk_id>
    # exemple: ref_and_rep_assemblies.acc_list_chunk.00
    file_name=$(basename -- "$chunk_file")
    db_name="${file_name%.*.*}"
    chunk_id="${chunk_file##*.}"

    formated_seq_fl=$out_db/${db_name}.chunk${chunk_id}.fna

    command=$command"module load devel/python/Python-3.11.1; python merge_and_format_assembly.py --assembly_selection assembly_summary_with_header.txt --assembly_root_dir /work/bank2/users_banks/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/ -o $formated_seq_fl $chunk_file -v;"

    command=$command'module load devel/Miniconda/Miniconda3; module load bioinfo/OBITools/1.2.11;'
    command=$command"obiconvert --fasta $formated_seq_fl --ecopcrdb-output=$out_db/${db_name}.chunk${chunk_id} -t ncbi_taxdump_2024_05 --nuc;"
    command=$command" gzip $formated_seq_fl"
    echo $command >> build_assemblies_ecopcr_db.sarray
done

# Launch
sarray --mem=200G build_assemblies_ecopcr_db.sarray

# Creation of name_seq with and without taxonomy. This allows me to have the species content in all my genomes and to keep this information for future comparison. 
## Retrieving the genome name using the taxid. This corresponds to the name of the sequences in the fastas.
zcat *fna.gz | grep ">" | cut -d' ' -f1-2 > name_seq_without_taxo.txt

## Use of Taxon kit, which takes the taxid as input and outputs the taxonomy down to species level.
cut -d'=' -f2 name_seq_without_taxo.txt | cut -d';' -f1 | taxonkit lineage –data-dir ../rpob/complete_and_chromosome/ncbi_tax_dumb/2024-15-04/ | taxonkit reformat –data-dir ../rpob/complete_and_chromosome/ncbi_tax_dumb/2024-15-04/ -f "{k};{p};{c};{o};{f};{g};{s}" -r "Unassigned" -P | cut -f 1,3-10 > name_seq.taxo

## Delete all those with not found IDs, otherwise obiconvert won't take them. You also need to remove the spaces in the names.

awk '{if (NF > 1 && $2 != "") print $0}' name_seq.taxo | sed 's/ /_/g' > filtered_name_seq.taxo

## Launch script to add the taxo to the sequence name
python complete_taxonomy.py name_seq_without_taxo.txt filtered_name_seq.taxo name_seq_with_taxo.txt

