# TaxonMarker

[TOC]

## Presentation
TaxonMarker is a tool developed in Python and Bash, designed to identify genetic markers specific to a given taxonomic rank. In metabarcoding sequencing analysis, the 16S rRNA gene is commonly used, although it has limitations, notably its presence in multiple copies and its low resolution at the level of closely related organisms. The aim of TaxonMarker is to select target genes that are potentially more effective than 16S at discriminating between species. The input is simply the desired taxonomic identifier (for example, Actinomycetes, taxid: 1760).

## Installation 

```bash=
git clone https://github.com/GTG1988A/TaxonMarker.git
```

## 1. formatting file bacterian OG

This step must only be performed by the programmer, each time OrthoDB is updated. It allows only bacterial ortholog groups to be retrieved from the OrthoDB database.

```bash=
cd TaxonMarker/1_formatting_file_bacterian_OG
python formatting_bacterial_orthologue_file_final.py -o ../Orthodb/odb11v0_OG2genes.tab -b only_line_bacteria.txt -u uniq_og_ids.txt -s ../Orthodb/odb11v0_level2species.tab -g ../Orthodb/odb11v0_OGs.tab -f Bacterial_OG.tab
```

## 2. search taxid and monocopy calculation

Retrieve the OGs containing the selected taxonomic rank. The identifiers of the species from this taxonomic rank are extracted, then the OGs containing at least one of these species are retained. This list of species is called identifiers_with_searchID_in_taxonomy.

This identifiers_with_searchID_in_taxonomy list is also used to check how many taxids in the ‘taxid’ column belong to this list. This allows us to determine the number of species belonging to our desired rank in the MT, as well as the percentage in relation to the total number of species.

For each MT corresponding to the specified identifiers, the script calculates and includes in the output tab :

    'OG_ID' : Ortholog group ID
    'ProteinCount' : Total number of proteins in the OG
    'SpeciesCount' : Number of unique species in the OG
    'nb_single_copy' : Number of proteins present in a single copy in at least one species
    'percent_single_copy' : Percentage of single-copy proteins in the OG
    'ProteinID' : List of all protein IDs. Protein IDs are made up of a number followed by ‘:’ followed by the taxid.
    'taxids' : List of taxonomic IDs corresponding to species in the OG
    'species' : List of species names corresponding to taxonomic IDs
    'TargetSpecies_Count' : Number of target species in the OG
    'TargetSpecies_Percentage' : Percentage of target species in the OG

With the ‘OG_selection.sh’ script you can select the best OGs


*Example command:*

```bash!
python search_taxid_and_monocopy_and_percentage_calculation.py -i ../1_formatting_file_bacterian_OG/Bacterial_OG.tab -f /BD_TaxonMarker/Orthodb/odb11v0_species.tab -l /BD_TaxonMarker/Orthodb/odb11v0_level2species.tab -s 1578 -o OG_1578.tab
```
### How OG_selection.sh works and Tips
The script can sort according to three criteria:

> -p <percent_single_copy> (mandatory, set to 0 if no selection is required) 
> [-c <TargetSpecies_Count>] Minimum number of species (because an OG containing only two species at 100% of our rank is less interesting than one containing 100) 
> [-t <TargetSpecies_Percentage>] Minimum percentage to include only OGs potentially corresponding to our rank of interest

We can find out how many species have descended from our taxonomic rank by using the following command (it is indicated in the --help of the script).
```bash=
$ grep -w "1578" /BD_TaxonMarker/Orthodb/odb11v0_levels.tab
1578    Lactobacillus   542691  14672   265
```

So we can try to go all out on the parameters 

```bash=
./OG_selection.sh -p  100 -c 265 -t 100 -o OG_1578_selected.tab OG_1578.tab
```

then gradually reduce until a satisfactory result is achieved.

## 3. fasta recovery

Downloading gene sequences in nucleic acid format. To do this, we use two APIs. This script retrieves the OGs selected in step 2. It extracts the protein ID of each protein in the OG, then uses the OrthoDB API to obtain the EMBL ID of the CDS. This ID is then used to download the nucleic sequence in FASTA format via the EMBL API. If no ID is found, this will be indicated in the logs.

*Example command:*
```bash!
python fastas_recovery.py ../2_search_taxid_and_monocopy_calculation/OG_1578_selected.tab
```

The script also adds the number of sequences contained in the OG fasta to the table.

## 4. primer design

This step contains several of them:

### a. Degeprime installation

```bash
git clone https://github.com/EnvGen/DEGEPRIME.git
```
### b. Alignment
Start sequence alignment of all OGs with Clustal Omega
*Example command:*
```bash!
mkdir -p alignment/

for files in ../3_fasta_recovery/*.fa; do
    file_name=$(basename "$fichier")
    clustalo -i "$files" -o "alignment/$file_name" --threads=4
done
```
### c. Degeprime trim launch

*NB: See more about trim in the degeprime github.*
*Example command:*
```bash!
for files in alignment/*.fa; do
    file_name=$(basename "$fichier" .fa)
    perl DEGEPRIME/TrimAlignment.pl -i "$files" -min 0.9 -o "alignment/trimmed_$file_name.fna"
done

```
### c. Génération de toutes les amorces possibles

Generation of all possible primers

All possible primers are generated by specifying the minimum and maximum primer length parameters, as well as the degeneracy values.

The command below generates a sarray file to run all possible parameter combinations.


*Example command:*


```bash!
#!/bin/bash

s_array_file='degeprime_multiple_params.sarray'

fasta_files='alignment/*.fa'  
result_dir='degeprime_result/'
aln_dir='alignment/'

degeneracies='8 12 24 48 96'
lengths='14 15 16 17 18 19 20 21 22 23 24'

echo '#' $fasta_files > $s_array_file
echo '# degeneracies' $degeneracies >> $s_array_file
echo '# lengths' $lengths >> $s_array_file
echo '# result:' $result_dir >> $s_array_file

for fasta_file in $fasta_files
do
  cog=$(basename "$fasta_file" .fa)
  trimmed_alignment="$aln_dir/trimmed_${cog}.fna"

  for d in $degeneracies
  do
    for l in $lengths
    do
      output_fl=$result_dir/${cog}_d${d}_l${l}.tsv

      echo "perl DEGEPRIME/DegePrime.pl -i $trimmed_alignment -d $d -l $l -o $output_fl" >> $s_array_file
    done
  done
done

echo sarray -J degeprime -o slurm_array_out/%j_%x.out $s_array_file
```


NB: Degeneracy is the total number of possible combinations of nucleotides that a degenerate primer can form. For example, for the ATCS primer, where S represents G or C, the possible combinations are ATCG and ATCC. The degeneracy is therefore equal to 2. Possible degeneracy values are 2, 4, 8, 16, 32, 64, 96, etc. Higher values are possible, but a degeneracy of 96 is already extremely unstringent.

### d. Concaténation des résultats

We concatenate the results to obtain one file per OGs

```bash!
#!/bin/bash

# Define the path of the directory containing the files
directory="degeprime_result"

# Recover file prefixes
prefixes=$(ls -1 "$directory" | grep -E ".*\.tsv$" | cut -d'_' -f1 | sort | uniq)

# Loop on each prefix

for prefix in $prefixes; do
    # Concatenate all files with the same prefix
    files=$(ls -1 "$directory"/"$prefix"* 2>/dev/null | grep -E ".*\.tsv$")
    if [ ! -z "$files" ]; then
        concatenated_file="$directory/concatenated_$prefix.tsv"
        cat $(echo "$files" | head -n1) > "$concatenated_file"
        for file in $(echo "$files" | tail -n+2); do
            tail -n +2 "$file" >> "$concatenated_file"
        done
    fi
done
```
### e. Adding metrics to primers

We run a program to add metrics to the primers. 

The statistics added are: minimum and maximum temperature, percentage of GC, presence of a GC clamp, percentage of number matching, whether the primer ends in a T, whether it is complementary to itself.

We will specify the desired maximum and minimum temperatures, as well as the minimum matching percentage (= number of sequences where a primer is found identically). Primers that do not meet these criteria will be deleted.

We will obtain a TSV file containing the primers matching our criteria, with the new statistics added to the table.

*Example command:*
```bash!
python process_primers_stat.py -i degeprime_result/concatenated_* -og ../3_fasta_recovery/updated_OG_selected_1578.tab -o result_stat_primers -nm 80 -tm_max 65 -tm_min 54
```

### f. Création du tableau des couples de primers.

Creation of a table of pairs of primers.

We run a program to generate the possible pairs of primers. We specify the minimum and maximum amplicon size between the two primers. Pairs that do not meet these criteria will be deleted. In addition, pairs with a temperature difference of more than 5 degrees between the primer forward and the primer reverse will also be eliminated.

*Example command:*
```bash!
#!/bin/bash
for file in result_stat_primers/concatenated_*; do echo "module load devel/python/Python-3.11.1;python couple_primer.py -i $file -f alignment/ --amplicon_min_size 150 --amplicon_max_size 590" >> tab_couple.sarray;done

sarray --meme=200G tab_couple.sarray
```

Each pair of primers will be given a total score, calculated by adding the scores for number matching and amplicon size. These parameters were deemed to be the most important after selection based on temperature and amplicon size.


### g. Concaténation des résultats et séléction des meilleurs couples.

Concatenate the results and select the best pairs.

We concatenate the results using concatenate_sort_result.sh and select the pairs we are interested in.

We will select the pairs with the three highest scores. We can then try to select those with the lowest degeneracy, a GC percentage of 50% or a GC clamp at the end, etc. However, it is important to note that the perfect primer does not necessarily exist.

A script called script_js_html_two_tab.py takes the list of selected primers as input. It generates a histogram showing the distribution of degeneracy. It is preferable for degeneracy to be more concentrated at the beginning of the primer, rather than at the end. This histogram therefore provides a better visualisation of this distribution.

TODO: graph of GC content? Position?


NB: If you run this program on a calculation cluster, you can use the primer_pipeline.sh script, adjusting the parameters inside. This allows you to run all the steps in one go.
## 5. EcoPCR

### a. Lancement de EcoPCR
Launching EcoPCR
EcoPCR is launched on a database formatted for the tool (see the section ‘The different TaxonMarker databases’).

This allows us to test our primers in silico. We can check whether they correctly capture the desired targets and whether they also capture other unexpected sequences.

EcoPCRs are usually run one at a time with the selected primers. (To be discussed if we will set up an automated process to run multiple EcoPCRs simultaneously).

To run an EcoPCR, use the script 1_commande_ecopcr.sh.
```bash=
./1_ecopcr_command.sh
```

It looks like this:

```bash=
#!/bin/bash
s_array_file='ecopcr_commands.sarray'

ecopcr_assembly_dir='/PATH/BD_TaxonMarker/BD_ecoPCR/ecoPCR_db/assemblies/*.adx'

ecopcr_outdir='primers'
mkdir -p $ecopcr_outdir
COG='58151at1578'
amorce_f='AAYATGGGKGTBGGNGAYG'
amorce_r='TTCWGGGAARAYBARYTGT'


for assembly in $ecopcr_assembly_dir
    do
    name=$(basename "$assembly" | sed 's/\.adx$//')

    echo "module load bioinfo/ecoPCR/1.0.1;ecoPCR -d $(echo "$assembly" | sed 's/\.adx$//')  $amorce_f $amorce_r  > primers/${name}.ecopcr;" >> $s_array_file
    done

```

You can specify a name for the output file by defining s_array_file. Once the first run has been completed, further analyses can be run.

The following information is also required:

    The path to the EcoPCR library (supplied by us)
    The name of the target gene or its OG
    The forward and reverse primers selected
Le script génère alors la commande EcoPCR que vous pouvez exécuter comme suit :

```bash=
sarray --mem=200G ecopcr_commands.sarray
```

The script then generates the EcoPCR command, which you can run as follows:

### b.  Analysis of EcoPCR results

#### 1_launch_format_ecopcr_result.sh

This script converts the EcoPCR results into FASTA format, adds taxonomy information and then launches the vsearch tool for sequence dereplication. It also removes all sequences containing amino acids, which is necessary for SWARM to work properly.

*Example command:*

```bash=
python /path/Lactobacillus/TaxonMarker/script_treatment_ecopcr_result/format_ecopcr_result.py -o . -t /path/BD_TaxonMarker/BD_ecoPCR/name_seq_with_taxo.txt primers/*ecopcr

```
It takes as an argument the name_seq_with_taxo.txt file in the EcoPCR_DB database for TaxonMarker (see the ‘Database’ section). This file contains the names of the sequences in FASTA format together with taxonomy information.

The final output file is named all_modified.fna.


#### 2_launch_swarm.sh

We use clustering with the Swarm tool to measure species discrimination. Swarm is a highly discriminating, robust and rapid tool, which groups amplicons with a high degree of similarity between them. We assume that two amplicons grouped in the same cluster come from the same species.

*Example command:*
```bash!
python /PATH/TaxonMarker/script_treatment_ecopcr_result/Launch_swarm.py -f all_modified.fna -s fichier_swarm.txt -o cluster.txt  -t 4 -a 1 -d 1
```
The output will be the file that swarm produces (ex:fichier_swarm.txt) and the file that formats swarm's results (ex:cluster.txt).

the format file provides all the taxonomic information for each cluster. 

#### 3_launch_stat_swarm.sh

This script allows us to obtain species discrimination measures specifically for our taxonomic rank. It offers the possibility of including or excluding certain terms in the taxonomy, in order to generate a clean statistical file, limited solely to our rank of interest.

##### Astuce au travers d'un cas concret:
Here's a practical tip:
Let's take the example of Lactobacillus:

All the members of this group no longer belong to the same genus as before, when they were all classified under the genus Lactobacillus. Since the new classification, several genera have been created from this genus. However, they still share the same family, Lactobacillaceae.


If I run my script including only the Lactobacillaceae family, i.e. using
-i f__Lactobacillaceae

I will have the file cluster_corrected.txt which only contains species belonging to the family Lactobacillaceae in their taxonomy, and the file rejected_clusters.txt with those that do not belong to this family.

However, when we look at the results, we see that the Lactobacillaceae family includes not only Lactobacillus, but also other genera that we want to exclude. We will then draw up a list of all the genera or species to be removed:

g__Leuconostoc g__Pediococcus g__Oenococcus g__Unassigned sp. g__Weissella g__Dellaglioa g__Periweissella g__Convivina


So we run our script again with the exclusion of these words

```bash=
python stats_report_taxo2.py -c cluster.txt  -o stats.txt -i f__Lactobacillaceae -e g__Leuconostoc g__Pediococcus g__Oenococcus g__Unassigned sp. g__Weissella g__Dellaglioa g__Periweissella g__Convivina -l rejected_clusters.txt -r cluster_corrected.txt

```

We will then obtain the file cluster_corrected.txt containing only species of the genus Lactobacillus. Statistics will be calculated from this file.

The rejected_clusters.txt file will contain all the species amplified by our primers that are not Lactobacillus.

Note: We are not obliged to apply specific inclusions and exclusions. We can choose to keep all results unchanged.

NB: For example, for the ‘Actinomycete’ class there will be no need for exclusion words. 


##### Analysing the results of the stats.txt file

This is how the script works:

Note: We are working on the complete taxonomy. If we only use the taxid, the results may vary, as two different taxids may have the same taxonomy if they correspond to strain taxids rather than species taxids. The script, however, uses a taxonomy that goes down to species level.

The script assesses whether a taxonomy is ‘well discriminated’ or ‘bad discriminated’ by examining the diversity of taxonomies present in each cluster.

For each cluster, the script extracts all the taxonomies associated with the sequences. If all the sequences in a cluster belong to the same taxonomy (i.e. there is only one unique taxonomy in the cluster), the cluster is considered to have ‘good discrimination’. The taxonomy is then added to a set of ‘well discriminated taxonomies’.

If a cluster contains several different taxonomies (i.e. several taxonomies are present in the same cluster), the cluster is considered to have ‘bad discrimination’. The script then adds the taxonomies present in this cluster to a set of ‘poorly discriminated taxonomies’.

The script uses sets to store good_discriminated_taxonomies and bad_discriminated_taxonomies.
Sets are data structures that do not allow duplicate elements to be stored. If a taxonomy is added several times to a set, it will only be counted once.

Ambiguous taxonomies: If a taxonomy appears both in well-discriminated clusters and in poorly-discriminated clusters, it is considered to be ‘ambiguous’ or to have mixed discriminations.



```bash=
Number of unique taxonomies: 291 #Total number of different taxonomies in clusters
Number of good discriminated taxonomies: 263
Percentage of good discriminated taxonomies: 90.38% 
Number of taxonomies both good and bad discriminated: 9
Total number of filtered clusters: 324 
Number of clusters with good discrimination: 304
Number of clusters with bad discrimination: 20
Percentage of clusters with good discrimination: 93.83%

# Below, we have the poorly discriminated clusters and the species that show both good and poor discrimination.
```

## 6. genome selection

We therefore have information on the discrimination of our target genes for the selected taxonomic rank. But isn't the 16S gene better in the end? To find out, we need to make a comparison.

During EcoPCR, we test our primers on a complete microbial database (see the section on the database). The first step is to filter this database to retain only those species corresponding to our taxonomic rank of interest. We select the genomes whose taxonomy contains the relevant criteria (for example, using the inclusion and exclusion words).

This will allow us to search these genomes for those that have a 16S sequence. In this way, we will be able to correctly compare our analysis of the target gene with that of 16S, since we will have selected the genomes of interest from the same database.

*Example command:*
For example, for Lactobacillus, we would run the search on our microbial database using the following criteria:

```bash=
python filter_seq_with_keywords.py -s /home/gagoutin/work/BD_TaxonMarker/BD_ecoPCR/name_seq_with_taxo.txt -o filtered_seq_with_taxo.txt -i f__Lactobacillaceae -e g__Leuconostoc g__Pediococcus g__Oenococcus g__Unassigned sp. g__Weissella g__Dellaglioa g__Periweissella g__Convivina
```


TODO: This script will output the total number of species in our database in the log. This information is interesting because we can compare it with the ‘Number of unique taxonomies:’ in our stats.txt file obtained just before. This allows you to see how many we catch in relation to the total. 

NB for lacto:
```bash=
gagoutin@genobioinfo1 ~/work/Lactobacillus/TaxonMarker/6_genome_selection $ cut -d" " -f3 filtered_seq_with_taxo.txt | sort -u | wc -l
365
```

he script will also produce the file genome_name_selected.txt, which contains a list of the names of the selected genomes. This file will allow us to filter the genomes of interest in our bank of 16S sequences, in order to retain only the 16S sequences corresponding to the desired ranks (e.g. Lactobacillus).


## 7. 16S comparison
We start by running the 16S_sequence_selected_extractor.py script, which takes our previously obtained list of genomes as input. This script extracts only the 16S sequences corresponding to the selected taxonomic rank.

```bash=
python 16S_sequence_selected_extractor.py --genome_file ../6_genome_selection/genome_name_selected.txt --tsv_file /PATH/BD_TaxonMarker/BD_16S/16S_with_taxonomy.tsv --fasta_file /PATH/BD_TaxonMarker/BD_16S/16S.fna --output_file filtered_16S.fna
```

We will have several outputs:

- The filtered fasta file: filtered_16S.fna
- name_seq_16S.txt: information on the taxonomy of each sequence
- process.log where at the end we will have the ‘Number of unique taxonomies’. It tells us the number of species that have a 16S

example: we have 365 different taxonomies for lactobacillus and we obtain 353 species with a 16S.


Format the fasta obtained with only the 16S sequences of the row selected with obiconvert.

```bash=
obiconvert --fasta filtered_16S.fna --ecopcrdb-output=ecoPCR_db/16S -t /PATH/BD_TaxonMarker/ncbi_tax_dumb/2024-15-04

```

Then run ecopcr
NB: Use the 16S primers commonly used in the laboratory. 
```bash=
ecoPCR -d ecoPCR_db/16S CCTACGGGNGGCWGCAG GACTACHVGGGTATCTAATCC > 16S.ecopcr
```


### Processing the results
The results are then processed as above:

#### 1_launch_format_ecopcr_result_16S.sh

```bash=
python script_treatment_ecopcr_result/format_ecopcr_result.py -o . -t name_seq_16S.txt 16S.ecopcr
```

#### 2_launch_swarm.sh
```bash!
python /home/gagoutin/work/Lactobacillus/TaxonMarker/script_treatment_ecopcr_result/Launch_swarm.py -f all_modified.fna -s fichier_swarm.txt -o cluster.txt  -t 4 -a 1 -d 1
```

### 3_launch_stat_swarm_16S.sh
```bash!
python stats_report_taxo2.py -c cluster.txt -t all_modified.fna -o stats.txt -l rejected_clusters.txt -r cluster_corrected.txt

```

NB: there's no need to include or exclude because we've filtered before.


#### Example of result analysis for Lactobacillus

```bash!
Number of unique taxonomies: 259
Number of good discriminated taxonomies: 159
Percentage of good discriminated taxonomies: 61.39%
Number of taxonomies both good and bad discriminated: 21
Total number of filtered clusters: 264
Number of clusters with good discrimination: 221
Number of clusters with bad discrimination: 43
Percentage of clusters with good discrimination: 83.71%

```

The clusters appear to be fairly well discriminated, but we note that many identical species are spread over several clusters. In fact, the well-discriminated taxonomies correspond to only 61%.


### Comparaison avec species_comparison_venn.py

```bash!
python species_comparison_venn.py --taxo_target_gene ../5_EcoPCR/uniq_taxo.txt --taxo_16S_gene uniq_taxo.txt --good_discrimination_target_gene ../5_EcoPCR/uniq_taxo_good_discriminated.txt --good_discrimination_16S uniq_taxo_good_discriminated.txt --output venn_diagram.png
```


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
/home/gagoutin/work/base_données//NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/GCA/018/630/415/GCA_018630415.1_ASM1863041v1_genomic.fna.gz GCA_018630415.1 1280
/home/gagoutin/work/base_données//NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/GCA/018/630/395/GCA_018630395.1_ASM1863039v1_genomic.fna.gz GCA_018630395.1 1280

# $ head genome_dirs.txt
/home/gagoutin/work/base_données/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/GCA/018/630/415
/home/gagoutin/work/base_données/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/GCA/018/630/395

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

#### 1. Récupérer l'assembly summary
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

*Utilisation du name_seq:*

Nous pouvons sélectionner les données qui nous intéressent en utilisant le même filtre que celui appliqué à la fin de TaxonMarker, à l'aide du script filter_seq_with_keywords.py. Par exemple, pour les Lactobacillus, il peut être difficile de les identifier tous car ils n'ont pas tous un caractère commun exclusif. Ils appartiennent tous à la famille Lactobacillaceae, mais il est possible que d'autres genres en fassent également partie.


Ce fichier peut également être utilisé pour ajouter des informations à un fichier FASTA. Si vous disposez uniquement des identifiants des séquences, vous pouvez retrouver les informations taxonomiques correspondantes grâce au fichier name_seq.txt. Ce fichier permet de faire le lien entre les identifiants des séquences et leurs informations taxonomiques, facilitant ainsi l'enrichissement de votre fichier FASTA avec les données de taxonomie pertinentes.

*Lancement de EcoPCR sur la base*

on lance ecoPCR comme ça: 

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
