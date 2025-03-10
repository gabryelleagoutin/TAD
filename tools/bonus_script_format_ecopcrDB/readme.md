# The different TaxonMarker databases
## OrthoDB

Database maintained by the University of Geneva and the SIB (Swiss Institute of Bioinformatics)

It contains groups of orthologous genes from the entire living kingdom
 = All the descendants of a particular single gene from the last common ancestor of these species.

This is the database used by the TaxonMarker tool to retrieve target genes. 

links to download OrthoDB:
```bash
wget -r https://data.orthodb.org/v11/download/
```


## 16S Database 
### Description de la base de données utilisées (identique pour EcoPCR)
Description of the database used (identical for EcoPCR)
We will extract the 16S from a database containing the NCBI microbial genomes (fasta and gff) contained in the RefSeq and GenBank databases with the GTDB taxonomy using the genome_updater tool.

This database comes from this download 
```bash
genome_updater.sh -g "archaea,bacteria" -d "refseq,genbank" -M "gtdb" -f "genomic.gff.gz,genomic.fna.gz" -a -o downloads -N -t 2 -V -Z -R 6 -L curl
```
with this tool:
https://github.com/pirovc/genome_updater
Bash script to download/update snapshots of files from NCBI genomes repository (refseq/genbank) with track of changes and without redundancy

It is kept at the Genotoul cluster in Toulouse, and all the information is available here:

```bash 
/work/bank2/users_banks/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/README
```

### 16S extraction

We need the scripts
- regions_analysis_fct.py
- identify_rrna_gene.py
- extract_16S23S.py

the extract_16S23S.py script needs two files
python extract_16S23S.py -a genome_dirs.txt -t genome_files_with_taxid.txt

genomes_dirs.txt contains the location of each genome fna and gff file
genome_files_with_taxid.txt contains the path to the genome file name, plus the name of the genome and its taxid 

```bash=
# $ head genome_files_with_taxid.txt
/PATH/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/GCA/018/630/415/GCA_018630415.1_ASM1863041v1_genomic.fna.gz GCA_018630415.1 1280
/PATH/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/GCA/018/630/395/GCA_018630395.1_ASM1863039v1_genomic.fna.gz GCA_018630395.1 1280

# $ head genome_dirs.txt
/PATH/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/GCA/018/630/415
/PATH/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/GCA/018/630/395

```

To get these files, you first need to run the scripts:
./1_find_dirs.sh
and 2_find_name.py, taking care to adjust the database paths used.

Once the extraction script has been run, we obtain two resulting files :

    16S.tsv: containing information about each 16S region extracted.
    16S.fna: containing the sequences of the 16S extracted.

I'm also going to add a 16S_with_taxonomy.tsv file which adds a taxonomy column to the 16S.tsv file. This file will allow us to select the taxonomies of interest when comparing the 16S with a target gene. For example, for Lactobacillus, we can filter to keep only these species and exclude the others.

To obtain these taxonomies, we use TaxonKit.

```bash=
cut -f 12 16S.tsv  | tail -n +2 | taxonkit lineage --data-dir ../ncbi_tax_dumb/2024-15-04/ | taxonkit reformat --data-dir ../ncbi_tax_dumb/2024-15-04/  -f "{k};{p};{c};{o};{f};{g};{s}" -r "Unassigned" -P | cut -f 1,3-10  > 16S.taxo

```

When the column is empty because there is no taxonomy, I add an Unassigned:
```bash=
awk 'BEGIN {FS=OFS="\t"} {if ($2 == "") $2 = "k__Unassigned;p__Unassigned;c__Unassigned;o__Unassigned;f__Unassigned;g__Unassigned;s__Unassigned"; print}' 16S.taxo > 16S_taxo_modified.txt

```

I remove unwanted spaces: 

```bash=
 sed -i 's/ /_/g' 16S_taxo_modified.txt
```

I keep only the second column and add the header ‘taxonomy’ to it

```bash=
cut -f2 16S_taxo_modified.txt | sed  '1s/^/taxonomy\n/' > col2.txt
```

I make sure that the files are encoded correctly:

```bash=
dos2unix col2.txt 16S.tsv
```
/!\ It's important to do this, because if the TSV file is encoded as ‘ASCII text, with CRLF line terminators’ instead of just ‘ASCII text’, the paste command won't work properly. You can check the encoding of the file using the following command:

```bash!
file 16S.tsv
```

and now I paste the taxonomy column into the tsv file:

```bash=
paste 16S.tsv col2.txt > 16S_with_taxonomy.tsv
```
## ncbi tax_dump

This gives us the ncbi taxonomy.
all the info can be found here:
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/

That's how you download it:
```bash=
wget -r https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/

```
## EcoPCR DataBase
### Description de la base de données utilisées (identique pour 16S)
Description of the database used (identical for 16S)
We are formatting a database containing the NCBI microbial genomes (fasta and gff) contained in the RefSeq and GenBank databases with the GTDB taxonomy using the genome_updater tool.

This database comes from this download:
```bash
genome_updater.sh -g "archaea,bacteria" -d "refseq,genbank" -M "gtdb" -f "genomic.gff.gz,genomic.fna.gz" -a -o downloads -N -t 2 -V -Z -R 6 -L curl
```
with this tool:
https://github.com/pirovc/genome_updater
Bash script to download/update snapshots of files from NCBI genomes repository (refseq/genbank) with track of changes and without redundancy

It is kept at the Genotoul cluster in Toulouse, and all the information is available here:
```bash 
/work/bank2/users_banks/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/README
```

That's how she's formatted:
### Formatting the database for EcoPCR

#### 1. Retrieving the summary assembly
```bash
cp  /work/bank2/users_banks/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/assembly_summary.txt .
```

#### 2.  Complete with the header

```bash
cat header.txt
assembly_accession      bioproject      biosample       wgs_master      refseq_category taxid   species_taxid   organism_name   infraspecific_name      isolate version_status  assembly_level  release_type genome_rep       seq_rel_date    asm_name        asm_submitter   gbrs_paired_asm paired_asm_comp ftp_path        excluded_from_refseq    relation_to_type_material       asm_not_live_date       assembly_typegroup    genome_size     genome_size_ungapped    gc_percent      replicon_count  scaffold_count  annotation_provider     annotation_name annotation_date total_gene_count        protein_coding_gene_count    non_coding_gene_count    pubmed_id

cat header.txt assembly_summary.txt> assembly_summary_with_header.txt
```

#### 3. separation into chunks to go faster because of the huge genome

```bash
nb_chunk=500
mkdir assembly_chunks
cut -f1 assembly_summary_with_header.txt > ref_and_rep_assemblies_acc_list.txt
split -a 3 -d --number=l/$nb_chunk ref_and_rep_assemblies_acc_list.txt assembly_chunks/ref_and_rep_assemblies.acc_list_chunk.
```


#### 4.the sarray command to format the database
```bash
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
```
#### 5. Launch sarray

```bash
sarray --mem=200G build_assemblies_ecopcr_db.sarray
```
#### 6. Creation of name_seq with and without taxonomy.. 
This allows me to have the species content in all my genomes and to keep this information for future comparison.

Retrieving the genome name using the taxid. This corresponds to the name of the sequences in the fastas.

```bash
zcat *fna.gz | grep ">" | cut -d' ' -f1-2 > name_seq_without_taxo.txt
```

Use of Taxon kit, which takes the taxid as input and outputs the taxonomy down to species level:
```bash
cut -d'=' -f2 name_seq_without_taxo.txt | cut -d';' -f1 | taxonkit lineage –data-dir ../rpob/complete_and_chromosome/ncbi_tax_dumb/2024-15-04/ | taxonkit reformat –data-dir ../rpob/complete_and_chromosome/ncbi_tax_dumb/2024-15-04/ -f "{k};{p};{c};{o};{f};{g};{s}" -r "Unassigned" -P | cut -f 1,3-10 > name_seq.taxo
```

 Delete all those with not found IDs, otherwise obiconvert won't take them. You also need to remove the spaces in the names

```bash
awk '{if (NF > 1 && $2 != "") print $0}' name_seq.taxo | sed 's/ /_/g' > filtered_name_seq.taxo
```

Launch script to add the taxo to the sequence name
```bash
python complete_taxonomy.py name_seq_without_taxo.txt filtered_name_seq.taxo name_seq_with_taxo.txt
```

#### 7. Utilisation 

*Use of name_seq:*

We can select the data we are interested in using the same filter as that applied at the end of TaxonMarker, using the script filter_seq_with_keywords.py. For example, in the case of Lactobacillus, it can be difficult to identify them all because they do not all have an exclusive common trait. They all belong to the Lactobacillaceae family, but it is possible that other genera are also included.


This file can also be used to add information to a FASTA file. If you only have the sequence identifiers, you can find the corresponding taxonomic information using the name_seq.txt file. This file provides a link between the sequence identifiers and their taxonomic information, making it easier to enrich your FASTA file with the relevant taxonomy data.

*Launch EcoPCR on the database*.

This is how you launch ecoPCR: 

```bash=
#!/bin/bash
s_array_file='ecopcr_commands.sarray'

ecopcr_assembly_dir='/PATH/BD_ecoPCR/ecoPCR_db/assemblies/*.adx'

ecopcr_outdir='primers'
COG='58151at1578'
amorce_f='AAYATGGGKGTBGGNGAYG'
amorce_r='TTCWGGGAARAYBARYTGT'


for assembly in $ecopcr_assembly_dir
    do
    name=$(basename "$assembly" | sed 's/\.adx$//')

    echo "module load bioinfo/ecoPCR/1.0.1;ecoPCR -d $(echo "$assembly" | sed 's/\.adx$//')  $amorce_f $amorce_r  > primers/${name}.ecopcr;" >> $s_array_file
    done

```
