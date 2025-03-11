# STEP 1: Gene selection
**conda: TaxonMarker_main environment**

You need to retrieve the file Bacterial_OG.tab, which contains only bacterial ortholog groups from the OrthoDB database and is formatted with one OG per line.

download link ? [TODO]
## 1.  search taxid and monocopy calculation

```bash=
cd 1_search_taxid_and_monocopy_calculation
```

Retrieve the OGs containing the selected taxonomic rank.
The identifiers of the species belonging to this taxonomic rank are extracted, then the OGs containing at least one of them are retained. This list of species is called identifiers_with_searchID_in_taxonomy.

This list is also used to check how many taxids in the ‘taxid’ column belong to this selection. This makes it possible to determine the number of species in the MT belonging to the desired rank, as well as their percentage of the total number of species.

For each OG corresponding to the specified identifiers, the script calculates and includes in the output tab:

    'OG_ID': Orthologous group ID
    'ProteinCount': Total number of proteins in the OG
    'SpeciesCount': Number of unique species in the OG
    ‘nb_single_copy’ : Number of proteins present in a single copy in at least one species
    ‘percent_single_copy’ : Percentage of single-copy proteins in the OG
    'ProteinID': List of all protein IDs. Protein IDs are composed of a number followed by ‘:’ followed by the taxid.
    'taxids': List of taxonomic IDs corresponding to species in the OG.
    'species': List of species names corresponding to the taxonomic IDs.
    ‘TargetSpecies_Count’ : Number of target species in the OG
    'TargetSpecies_Percentage' : Percentage of target species in the OG

The ‘OG_selection.sh’ script can be used to select the best OGs.

Example:
```bash!
python search_taxid_and_monocopy_and_percentage_calculation.py -i ../../../data_test/0_File_bacterian_OG/test_Bacterial_OG_small.tab -f $PATH_ORTHODB/Orthodb/odb11v0_species.tab -l $PATH_ORTHODB/Orthodb/odb11v0_level2species.tab -s 1578 -o test_output_OG_1578_home.tab
```

NB :
- The data_test folder contains the results expected when you run the command with the test data. Feel free to check them.
- The launch_search_taxid_and_monocopy_and_percentage_calculation.sh script is designed to be used on a calculation cluster. You can adapt it to suit your needs.

### How OG_selection.sh works and Tips
The script can sort according to three criteria:

> -p <percent_single_copy> (mandatory, set to 0 if no selection is desired). 
> [-c <TargetSpecies_Count>] Minimum number of species (because an OG containing only two species at 100% of our rank is less interesting than one containing 100). 
-t <TargetSpecies_Percentage>] Minimum percentage to include only OGs potentially corresponding to our rank of interest.

We can already find out how many species are descending from our taxonomic rank by using the following command (it is indicated in the --help of the script).

```bash=
$ grep -w "1578" $PATH_ORTHODB/Orthodb/odb11v0_levels.tab
1578    Lactobacillus   542691  14672   265
# it's the last column, 265 species.
```

So we can try to be as stringent as possible:
```bash=
./OG_selection.sh -p  100 -c 265 -t 100 -o OG_1578_selected.tab OG_1578.tab
```

then gradually reduce until a satisfactory result is achieved.

Example:
```bash!
./OG_selection.sh -p 0 -c 250 -t 100 -o test_output_OG_1578_selected_home.tab test_output_OG_1578_home.tab
```

## 2. fasta recovery

Downloading gene sequences in nucleic acid format

To do this, we use two APIs. This script retrieves the OGs selected in step 2, extracts the protein ID of each protein in the OG, and then queries the OrthoDB API to obtain the EMBL ID of the CDS. This ID is then used to download the nucleic sequence in FASTA format via the EMBL API.

If no ID is found, this will be indicated in the logs.

Example:
```bash!
cd 2_fasta_recovery

python fastas_recovery.py ../1_search_taxid_and_monocopy_calculation/test_output_OG_1578_selected_home.tab
```

The script also adds the number of sequences contained in the OG FASTA file to the table.

The output is :

- The FASTA files for each OG
- The updated table (in our example: updated_test_output_OG_1578_selected_home.tab)
- An HTML file providing information on both the orthologous gene group (OG) and the gene itself

You have obtained your orthologous genes. Now you need to move on to step 2 


