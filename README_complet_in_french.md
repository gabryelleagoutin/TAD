# TaxonMarker 

## Présentation

TaxonMarker est un outil développé en Python et Bash, conçu pour identifier des marqueurs génétiques spécifiques à un rang taxonomique donné. Lors de l'analyse de séquençage métabarcoding, le gène de l'ARNr 16S est couramment utilisé, bien qu'il présente des limites, notamment sa présence en copies multiples et sa faible résolution au niveau des organismes étroitement apparentés. Le but de TaxonMarker est de selectionner des gènes cibles potentiellement plus efficaces que le 16S pour discriminer les espèces. En entrée, il suffit de fournir l'identifiant taxonomique souhaité (par exemple, famille "Actinomycetes", taxid: 1760).


**TaxonMarker est divisé en 3 pipelines:**

* ETAPE 1: Selection des gènes 
Au cours de cette étape, nous allons sélectionner, dans la base de données OrthoDB, des groupes de gènes orthologues qui incluent un nombre suffisant d'espèces correspondant au rang taxonomique défini. Nous obtiendrons un fasta par groupe orthologue. Nous obtenons aussi une fiche récapitulative de toute les informations que nous avons sur le gène. 

  Voir le dossier ... pour plus de détail et toute les étapes.

* ETAPE 2 : Design des primers 
Durant cette étape, nous alignons les séquences FASTA afin de générer des paires de primers. Nous les sélectionnons ensuite en fonction de divers critères, tels que la température de fusion, la présence d'un GC clamp, etc. Nous obtenons les meilleurs couples possible avec une sortie graphique de leurs caractéristiques.

   Voir le dossier ... pour plus de détail et toute les étapes.

* ETAPE 3 : Validation in silico des primers et des amplicons.
Nous réalisons une PCR in silico pour valider nos amorces, puis nous menons une analyse discriminatoire à l’échelle de l’espèce afin de vérifier si l’amplicon ciblé est suffisamment discriminant. Enfin, pour les espèces détectées par la PCR in silico, nous récupérons leurs séquences d’amplicons connues (par exemple le 16S) et appliquons le même traitement discriminatoire pour déterminer si notre nouveau couple de primers se révèle plus performant.

  Voir le dossier ... pour plus de détail et toute les étapes.

## Installation de TaxonMarker

```
git clone https://github.com/GTG1988A/TaxonMarker.git
cd TaxonMarker
```

## Préparation de l'environnement

### Création des différents environnements conda

```bash=
conda env create --name TaxonMarker_main --file environment_main.yml
conda activate TaxonMarker_main
```

Pour les étapes de PCR in silico et l'utilisation d'outils de la suite Obitools, une version de python 2.7 doit être utilisé. Il faut donc utiliser un second environnement dédié: 


```bash=
conda env create --name TaxonMarker_ecopcr --file environment_ecopcr.yml
conda activate TaxonMarker_ecopcr
```


### Téléchargement de la base de données OrthoDB 

```bash=
wget -r -np  --cut-dirs=2  https://data.orthodb.org/v12/download/
```

## Lancement de l'étape 1
Tout au long de ces étapes, des exemples de commandes sont illustrés à l’aide d’un jeu de données de test disponible dans le dossier ```data_test.```
=> aller dans le dossier STEP1_GENES_SELECTION
```
cd STEP1_GENES_SELECTION
```

# Descriptions des bases de données disponible pour TaxonMarker: 

**OrthoDB est la seule base de données requise pour la génération de vos primers.**

Cependant, pour les tests de PCR in silico ainsi que la comparaison avec le 16S, l'utilisation de différentes bases de données sera nécessaire.

Les étapes de traitement permettant de les convertir au format compatible avec TaxonMarker, ainsi que les scripts associés, se trouvent dans le dossier tools/bonus_script_format_ecopcrDB.


## OrthoDB

https://www.orthodb.org/

Base de données maintenue par l’Université de Genève et le SIB (Swiss institue of Bioinformatics)

Elle contient des groupes de gènes orthologues de tout le règne vivant
 = Tous les descendants d'un gène unique particulier du dernier ancêtre commun de ces espèces

## EcoPCR DataBase
Nous formatons une base de données contenant **tout** les Génomes microbiens NCBI (fasta et gff) contenus dans les bases de données RefSeq et GenBank avec la taxonomie GTDB via l'outil genome_updater (https://github.com/pirovc/genome_updater).

Fichier important: name_seq_with_taxo.txt. 
Cela me permet d'avoir le contenu de l'espèce dans tous mes génomes et de conserver cette information pour des comparaisons futures. Vous en aurez besoin pour le lancement de certain script.

Elle est conservé au cluster Genotoul de Toulouse.
Taille: 781 Go

## 16S Database 

Nous avons extrait les séquences d'ARNr 16S à partir d'une base de données contenant les Génomes microbiens NCBI (fasta et gff) contenus dans les bases de données RefSeq et GenBank avec la taxonomie GTDB avec l'outil genome_updater. 

*NOTE : C'est donc exactement le même téléchargement que pour la banque complète d'EcoPCR.*


Elle n'est pas au format accepté par EcoPCR. Vous devez formater votre selection de séquences. Tout est indiqué dans 

```tools/STEP3_IN_SILICO_VALIDATION/2_Amplicon_comparison/```

Fichier important: 16S_with_taxonomy.tsv 
Cela me permet d'avoir le contenu de l'espèce dans tous mes génomes et de conserver cette information pour des comparaisons futures. Vous en aurez besoin pour le lancement de certain script.

Elle est conservé au cluster Genotoul de Toulouse.
Taille : 1.5 Go 
## ncbi tax_dump

Cela nous permets d'avoir la taxonomie du ncbi.
toute les infos se trouvent ici:
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/

On le télécharge comme ceci:
```bash=
wget -r https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/
```

Taille: 1.9 Go

______________________________________________________________________________

[]

# 0. File bacterian OG
Cette étape doit être effectuée uniquement par le programmeur à chaque mise à jour microbienne d'OrthoDB. Elle permet d'extraire uniquement les groupes orthologues bactériens de la base de données OrthoDB.

```bash=
python formatting_bacterial_orthologue_file_final.py -o $PATH_ORTHODB/Orthodb/odb11v0_OG2genes.tab -b only_line_bacteria.txt -u uniq_og_ids.txt -s ../Orthodb/odb11v0_level2species.tab -g ../Orthodb/odb11v0_OGs.tab -f Bacterial_OG.tab
```

------------------------------------------------------------------------------
[]
# ETAPE 1: Selection des gènes
*Environnement conda: TaxonMarker_main* 

Vous devez récupérer le fichier Bacterial_OG.tab. Il contient uniquement les groupes orthologues bactériens de la base de données OrthoDB et est formaté avec un OG par ligne.

lien de téléchargement ? [TODO]

## 1.  search taxid and monocopy calculation

```bash=
cd 1_search_taxid_and_monocopy_calculation
```

Récupération des OGs contenant le rang taxonomique sélectionné.
On extrait les identifiants des espèces appartenant à ce rang taxonomique, puis on conserve les OGs qui en contiennent au moins une. Cette liste d'espèces est appelée identifiers_with_searchID_in_taxonomy.

Cette liste est également utilisée pour vérifier combien de taxid dans la colonne "taxid" appartiennent à cette sélection. Cela permet de déterminer le nombre d'espèces de l'OG appartenant au rang souhaité, ainsi que leur pourcentage par rapport au nombre total d'espèces.

Pour chaque OG correspondant aux identifiants spécifiés, le script calcule et inclut dans l'onglet de sortie :

    'OG_ID' : ID du groupe orthologue
    'ProteinCount' : Nombre total de protéines dans l'OG
    'SpeciesCount' : Nombre d'espèces uniques dans l'OG
    'nb_single_copy' : Nombre de protéines présentes en copie unique dans au moins une espèce
    'percent_single_copy' : Pourcentage de protéines en copie unique dans l'OG
    'ProteinID' : Liste de tout les ID protéiques. Les ID des protéines sont composé de chiffre puis de ":" puis du taxid.
    'taxids' : Liste des IDs taxonomiques correspondant aux espèces dans l'OG
    'species' : Liste des noms d'espèces correspondant aux IDs taxonomiques
    'TargetSpecies_Count' : Nombre d'espèces cibles dans l'OG
    'TargetSpecies_Percentage' : Pourcentage d'espèces cibles dans l'OG

Avec le script "OG_selection.sh" on peut selectionner les meilleurs OGs


*Exemple de commande:*

```bash!
python search_taxid_and_monocopy_and_percentage_calculation.py -i ../../../data_test/0_File_bacterian_OG/test_Bacterial_OG_small.tab -f $PATH_ORTHODB/Orthodb/odb11v0_species.tab -l $PATH_ORTHODB/Orthodb/odb11v0_level2species.tab -s 1578 -o test_output_OG_1578_home.tab
```

NB :
- Le dossier data_test contient les résultats attendus lorsque vous exécutez la commande avec les données de test. N'hésitez pas à les vérifier.
- Le script launch_search_taxid_and_monocopy_and_percentage_calculation.sh est conçu pour être utilisé sur un cluster de calcul. Vous pouvez l'adapter en fonction de vos besoins.

### fonctionnement de OG_selection.sh et Astuces
Le script peut trier selon trois critères :

> -p <percent_single_copy> (obligatoire, mettre à 0 si aucune sélection n'est souhaitée) 
> [-c <TargetSpecies_Count>] Nombre minimum d'espèces (car un OG contenant seulement deux espèces à 100 % de notre rang est moins intéressant qu'un OG en contenant 100) 
> [-t <TargetSpecies_Percentage>] Pourcentage minimum pour n'inclure que les OGs correspondant potentiellement à notre rang d'intérêt

On peut déjà savoir combien d'espèces descendent de notre rang taxonomique en utilisant la commande suivante (elle est indiquée dans le --help du script).

```bash=
$ grep -w "1578" $PATH_ORTHODB/Orthodb/odb11v0_levels.tab
1578    Lactobacillus   542691  14672   265
# C'est la dernière colonne : 265 espèces.
```

On peut donc essayer d'être le plus stringeant possible :

```bash=
./OG_selection.sh -p  100 -c 265 -t 100 -o OG_1578_selected.tab OG_1578.tab
```

puis réduire progressivement jusqu'à obtenir un résultat satisfaisant.

exemple:
```bash!
./OG_selection.sh -p 0 -c 250 -t 100 -o test_output_OG_1578_selected_home.tab test_output_OG_1578_home.tab
```

## 2. fasta recovery

Téléchargement des séquences des gènes au format acide nucléique

Pour cela, nous utilisons deux API. Ce script récupère les OGs sélectionnés à l'étape 2, extrait l'ID protéique de chaque protéine dans l'OG, puis interroge l'API d'OrthoDB pour obtenir l'ID EMBL du CDS. Cet identifiant est ensuite utilisé pour télécharger la séquence nucléique au format FASTA via l'API EMBL.

Si aucun ID n'est trouvé, cela sera indiqué dans les logs.


*Exemple de commande:*
```bash!
cd 2_fasta_recovery

python fastas_recovery.py ../1_search_taxid_and_monocopy_calculation/test_output_OG_1578_selected_home.tab
```

Le script ajoute également dans le tableau le nombre de séquences contenues dans le fichier FASTA de l'OG.

En sortie, vous obtenez :

- Les fichiers FASTA pour chaque OG
- Le tableau mis à jour (dans notre exemple : updated_test_output_OG_1578_selected_home.tab)
- Un fichier HTML fournissant des informations à la fois sur le groupe de gènes orthologues (OG) et sur le gène lui-même



Vous avez obtenu vos gènes orthologues. Il faut maintenant passer à l'étape 2 [ici]


_________________________________________________________________________________________________________________

# STEP2_PRIMER_DESIGN

## Installation de Degeprime 

Vous avez besoin d'installer cet outil. 

```bash
git clone https://github.com/EnvGen/DEGEPRIME.git
```

## Lancement de 1_primer_pipeline.sh
Pour les utilisateurs d'un cluster de calcul, ce script lance toute les étapes d'un coup. Pour cela, vous devez faire une commande comme ceci:

```bash!
sbatch -J primer_pipeline  -o %x.out -e %x.err -p workq --mem=80G -c 4 --wrap="sh 1_primer_pipeline.sh --fasta_path ../STEP1_GENES_SELECTION/2_fasta_recovery/ --degeneracies_values 72 96 --min_primer_length 14 --max_primer_length 24 --tab_og_updated ../STEP1_GENES_SELECTION/2_fasta_recovery/test_output_updated_test_OG_1578_selected.tab --number_matching_max 80 --tm_max 65 --tm_min 54 --amplicon_min_size 150 --amplicon_max_size 590"
```

descriptions des paramètres:
```bash
Options:
  --fasta_path <path>            Path to the directory containing the .fa files (mandatory, no default value).
  --value_degeneracies <int>     List of degeneracy values (e.g., 48 96). Mandatory.
  --min_primer_length <int>      Minimum primer length. Default: 14.
  --max_primer_length <int>      Maximum primer length. Default: 24.
  --tab_og_updated <path>        File for the -og parameter in process_primers_stat.py. This is the updated OG file after downloading the FASTA files. (default is empty, so it is mandatory if Step 6 is used).
  --number_matching_max <int>    Parameter -nm (default: 80).is the number of sequences that match this primer divided by the total number of sequences * 100. This is the minimum we expect.
  --tm_max <int>                 Parameter -tm_max (default: 65).
  --tm_min <int>                 Parameter -tm_min (default: 54).
  --amplicon_min_size <int>      Minimum amplicon size (default: 150).
  --amplicon_max_size <int>      Maximum amplicon size (default: 590).
  -h, --help                     Display this message and exit.
```

En sortie, vous obtenez un fichier sorted_results.tsv contenant les couples d'amorces avec les trois meilleurs scores.

Pour comprendre chaque étape ou si vous n'êtes pas sur un cluster de calcul, veuillez lire la suite.

## Descriptions des différentes étapes et lancement sans 1_primer_pipeline.sh

### a. Lancement de l'alignement des séquences de tout les OGs avec Clustal Omega
**conda environement TaxonMarker_main**

*Exemple de commande:*
```bash!
mkdir -p alignment/

for files in ../STEP1_GENES_SELECTION/2_fasta_recovery/*.fa; do
     file_name=$(basename "$files")
     clustalo -i "$files" -o "alignment/$file_name" --threads=4
done

```
### b. Lancement de trim de Degeprime 

*NB: Voir + à propos de trim dans le github de degeprime.*
*Exemple de commande:*
```bash!
for files in alignment/*.fa; do
    file_name=$(basename "$files" .fa)
    perl DEGEPRIME/TrimAlignment.pl -i "$files" -min 0.9 -o "alignment/trimmed_$file_name.fna"
done

```
### c. Génération de toutes les amorces possibles

Toutes les amorces possibles sont générées en spécifiant les paramètres de longueur minimale et maximale du primer, ainsi que les valeurs de  dégénérescences.


*Exemple de commande:*


```bash!
fasta_files='alignment/*.fa'  
result_dir='degeprime_result/'
aln_dir='alignment/'

degeneracies='72 96' #a ajuster ici 
lengths='17 18 19 20' #a ajuster ici 

mkdir -p "$result_dir"

for fasta_file in $fasta_files
do
  cog=$(basename "$fasta_file" .fa)
  trimmed_alignment="$aln_dir/trimmed_${cog}.fna"

  for d in $degeneracies
  do
    for l in $lengths
    do
      output_fl=$result_dir/${cog}_d${d}_l${l}.tsv

      echo "Launch : perl DEGEPRIME/DegePrime.pl -i $trimmed_alignment -d $d -l $l -o $output_fl"
      perl DEGEPRIME/DegePrime.pl -i "$trimmed_alignment" -d "$d" -l "$l" -o "$output_fl" > degeprime.log
    done
  done
done
```

NB: La dégénérescence correspond au nombre total de combinaisons possibles de nucléotides qu'un primer dégénéré peut former. Par exemple, pour le primer ATCS, où S représente G ou C, les combinaisons possibles sont ATCG et ATCC. La dégénérescence est donc égale à 2. Les valeurs possibles de dégénérescence peuvent être 1, 2, 3, 4, 6, 8, 9, 12, 16, 18, 24, 27, 32, 36, 48, 54, 64, 72, 81, 96, 108 etc. 


```
-d should be a possible degeneracy, i.e. 1, 2, 3, 4, 6, 8, 9, 12, and so forth (or more generally, a number > 0 that can be expressed as 2^i * 3^j, where i and j are integers or 0).
```


### d. Concaténation des résultats

On concatène les résultats pour avoir un fichier par OGs

```bash!
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

### e. Ajout de métrique sur les primers

Nous lançons un programme qui ajoutera des métriques sur les primers.

Les statistiques ajoutées sont: température minimale et maximale, pourcentage de GC, présence d'un GC clamp, pourcentage de number matching, si le primer fini par un T, si il est complémentaire à lui même.

Nous spécifierons les températures maximales et minimales souhaitées, ainsi que le pourcentage minimum de correspondance (= nombre de séquences où une amorce est trouvée de manière identique). Les primers ne répondant pas à ces critères seront supprimés.

Nous obtiendrons un fichier TSV contenant les primers correspondant à nos critères, avec les nouvelles statistiques ajoutées dans le tableau.

*Exemple de commande:*

```bash!
mkdir result_stat_primers
python process_primers_stat.py -i degeprime_result/concatenated_* -og ../STEP1_GENES_SELECTION/2_fasta_recovery/updated_test_output_OG_1578_selected_home.tab -o result_stat_primers -nm 80 -tm_max 65 -tm_min 54
```

### f. Création du tableau des couples de primers.

Nous lançons un programme qui générera tout les couples de primers possibles. Nous spécifierons la taille minimale et maximale de l'amplicon entre les deux primers. Les couples ne respectant pas ces critères seront supprimés. De plus, les couples avec un écart de température supérieur à 5 degrés entre le primer forward et le primer reverse seront également éliminés.

*Exemple de commande:*
```bash!
for file in result_stat_primers/concatenated_*; do
    echo "Treatment of : $file"
    python couple_primer.py -i "$file" -f alignment/ --amplicon_min_size 150 --amplicon_max_size 590
done

```

Chaque couple de primers se verra attribuer un score total, calculé en additionnant le score pour le number matching et celui pour la taille de l'amplicon. Nous avons jugé que ces paramètres étaient les plus importants après avoir effectué la sélection basée sur les températures et les tailles d'amplicon.


### g. Concaténation des résultats et séléction des meilleurs couples.

On concatène les résultats avec 0_script_step_by_step/concatenate_sort_result.sh et on obtient le fichier sorted_results.tsv contenant les couples ayant les trois couples meilleurs scores. 

Par la suite, nous pourrons essayer de choisir ceux ayant la dégénérescence la plus faible, un pourcentage de GC de 50 % ou un GC clamp à la fin, etc.


## Visualisation 
Un script nommé primer_metrics_visualization.py prend en entrée le fichier sorted_results.tsv et génère :

    Par couple de primers :
        Une fiche récapitulative du couple de primers avec toutes ses métriques.
        Un graphique montrant la répartition de la dégénérescence. En effet, il est préférable que la dégénérescence soit plus concentrée au début du primer plutôt qu'à la fin. Cet histogramme permet donc une meilleure visualisation de cette répartition.
        Un graphique montrant la concentration en GC aux différentes positions.

    Par OG :
        Un graphique de la position des primers sur l'alignement du gène.

Dans le cadre de analyse avec l'outils FROGS, nous avons également besoin des tailles théoriques minimale et maximale de l'amplicon. Pour cela, nous lançons des étapes de PCR in silico afin de récupérer cette information. Cette étape peut être effectuée avant d'exécuter primer_metrics_visualization.py.

Toutes les étapes, si vous êtes sur un cluster de calcul, sont automatisées via les scripts suivants :

    1_Obiconvert_fasta_OG.sh
    2_launch_ecopcr_and_add_amplicon_length_info.sh, qui exécute ecopcr_and_add_amplicon_length_info.py
    3_launch_primer_metrics_visualization.sh, qui exécute primer_metrics_visualization.py

Vous n'avez qu'à adapter les fichiers .sh selon vos besoins.

**Sinon, voici les différentes étapes:**


## a. Obiconvert
***Il faut activé l'environnement conda TaxonMarker_ecopcr***

Ce script permet de formater les fichiers FASTA dans un format compatible avec ecoPCR.

Vous aurez besoin de la banque taxdump. Pour en savoir plus, consultez le README d'accueil, dans la section "Description des bases de données disponibles pour TaxonMarker".


```bash!
Usage: ./1_Obiconvert_fasta_OG.sh /path/to/ncbi_tax_dump fasta_file1 [fasta_file2 ...] or all fasta *.fa.
  the best use is to select the OG fastas in the results table.
  Formatting is a fairly long process, so select the ones you have rather than formatting them all.

```

exemple de commande
```
./1_Obiconvert_fasta_OG.sh $PATH_NCBI_TAX_DUMP/ ../../STEP1_GENES_SELECTION/2_fasta_recovery/60099at1578_fasta.fa
```

# b.ecoPCR 
On exécute le script ecopcr_and_add_amplicon_length_info.py, qui prend en entrée le fichier sorted_results.tsv contenant tous les meilleurs couples de primers.

Il effectue une ecoPCR pour chaque couple, détermine les tailles minimale et maximale de l'amplicon, puis met à jour sorted_results.tsv avec ces informations.

exemple de commande
```
python 2_launch_ecopcr_and_add_amplicon_length_info.py ../sorted_results.tsv ecoPCR_db_*
``` 

# c. Lancement de la visualisation
***Environnement conda TaxonMarker_main****

/!\ il faut activer l'environnement conda TaxonMarker_main. Sinon, vous aurez des erreurs  car ce n'est pas la bonne version de python.


example de commande

```
python 3_primer_metrics_visualization.py -i ../sorted_results_updated.tsv -o primer_metrics_visualization.html
```


Il faut maintenant sélectionner les couples d'amorces qui vous semblent les plus appropriés. Vous avez terminé et obtenu vos amorces.

Si vous souhaitez les tester, rendez-vous sur la page d'accueil pour lire les informations dans la section "Description des bases de données disponibles pour TaxonMarker", puis passez à l'étape 3.

______________________________________________________________________

# STEP3_IN_SILICO_VALIDATION

## 1. Validation de nos amorces
Lancement d'ecoPCR sur une base de données formatée pour l'outil

Cette étape permet de tester les amorces in silico. Elle permet de vérifier si elles capturent correctement les cibles souhaitées et si elles amplifient également d'autres séquences non prévues.

### a. Lancement de EcoPCR
***Environnement conda TaxonMarker_ecopcr***


Pour lancer une EcoPCR, utilisez le script 1_commande_ecopcr.sh.

```bash=
./1_ecopcr_command.sh -h 

Usage: ./1_ecopcr_command.sh -a '<primer list>' -d '<path to .adx files>'
  -a   Comma-separated list of primer pairs, each consisting of a forward primer and a reverse primer. Primers are separated by a space and surrounded by an inverted comma. e.g. ‘TSRTCAAGAACRTBGARR SAYGTYCTGBACRTCRTC,TSRTCAAGAACRTBGARR AYGTYCTGNACRTCGTCV’.
  -d   Access path to the database for EcoPCR, the folder containing the .adx files. e.g. /BD_TaxonMarker/BD_ecoPCR/ecoPCR_db/assemblies/

```

Ce script génère alors les commandes EcoPCR. Pour le lancer:

**- Si vous êtes sur un cluster**

exemple de commande

```bash!
./1_ecopcr_command.sh -a 'GTNCCDCAYGGYGGYGG GCTTCNARDGCCCADACTT' -d $PATH_BD_ECOPCR_ASSEMBLIES/
```

Puis executer la commande comme ceci:

```bash=
sarray --mem=200G ecopcr_commands.sarray
```

**- Si vous n'êtes pas sur le cluster**

Commentez la ligne 59 du script 1_ecopcr_command.sh : 
```
# echo "module load devel/Miniconda/Miniconda3; source $MY_CONDA_PATH; conda activate TaxonMarker_EcoPCR; ecoPCR -d $(echo "$assembly" | sed 's/\.adx$//') $primer_f $primer_r > $ecopcr_outdir/$outdir/${name}.ecopcr;" >> $s_array_file
```

Et decommentez la ligne 61 du script 1_ecopcr_command.sh:
```
echo "ecoPCR -d $(echo "$assembly" | sed 's/\.adx$//') $primer_f $primer_r > $ecopcr_outdir/$outdir/${name}.ecopcr;" >> $s_array_file
```

Lancer 1_ecopcr_command.sh
```
./1_ecopcr_command.sh -a 'GTNCCDCAYGGYGGYGG GCTTCNARDGCCCADACTT' -d $PATH_BD_ECOPCR_ASSEMBLIES/
```

Activer l'éxecution de votre script 
```
 chmod +x ecopcr_commands.sarray
```
Activez l'environnement conda TaxonMarker_ecopcr puis lancer cette commande

```
./ecopcr_commands.sarray
```

### b. Analyse de la discrimination Taxonomique

***Environnement conda TaxonMarker_main***

#### 1_launch_format_ecopcr_result.sh

Ce script permet de convertir les résultats d'ecoPCR au format FASTA, d'ajouter les informations de taxonomie, puis de lancer l'outil vsearch pour la déréplication des séquences. Il supprime également toutes les séquences contenant des acides aminés, ce qui est nécessaire pour que SWARM fonctionne correctement.

Nous conservons les informations des séquences supprimées par la déréplication dans le fichier output_vsearch_cluster.txt. Cela nous permet de les ajouter ensuite aux résultats de SWARM, afin de ne pas caractériser une espèce comme correctement discriminée alors que son équivalent a été supprimé lors de la déréplication.

*Exemple de commande:*

```bash=
python ../script_treatment_ecopcr_result/format_ecopcr_result.py -o . -t $PATH_BD_ECOPCR/name_seq_with_taxo.txt result/1_GTNCCDCAYGGYGGYGG-GCTTCNARDGCCCADACTT/*ecopcr

```

Le fichier fasta de sorti se nomme **all_modified.fna**

#### 2_launch_swarm.sh

Nous utilisons le clustering avec l'outil Swarm pour mesurer la discrimination des espèces. Swarm est un outil hautement discriminant, robuste et rapide, qui regroupe les amplicons présentant de fortes similitudes.

Nous partons du principe que deux amplicons regroupés dans le même cluster proviennent de la même espèce.


*Exemple de commande:*
```bash!
python ../script_treatment_ecopcr_result/Launch_swarm.py -f all_modified.fna -s fichier_swarm.txt -o cluster.txt  -t 4 -a 1 -d 1 -vsearch output_vsearch_cluster.txt -taxo $PATH_BD_ECOPCR/name_seq_with_taxo.txt
```

En sortie, nous obtenons :

- Le fichier brut généré par Swarm (fichier_swarm.txt).
- Un fichier formaté contenant les résultats de Swarm (cluster.txt), rendant les données plus lisibles avec l'intégration de toutes les taxonomies.

#### taxonomic_discrimination_at_species.py

Ce script permet d'obtenir des mesures de discrimination des espèces, soit spécifiquement pour un rang taxonomique donné, soit de manière plus générale. Il offre la possibilité d'inclure ou d'exclure certains termes dans la taxonomie afin de générer un fichier statistique propre, limité uniquement au rang d'intérêt.


Voici les paramètres: 

```
  -h, --help            show this help message and exit
  -c CLUSTER_FILE, --cluster_file CLUSTER_FILE
                        Input file containing cluster information.
  -o OUTPUT_STATS_FILE, --output_stats_file OUTPUT_STATS_FILE
                        Output file for the statistics report.
  -r OUTPUT_CORRECTED_FILE, --output_corrected_file OUTPUT_CORRECTED_FILE
                        Output file for the filtered clusters.
  -i [INCLUDE_KEYWORDS ...], --include_keywords [INCLUDE_KEYWORDS ...]
                        Keywords that must be present in taxonomy to keep a line.
  -e [EXCLUDE_KEYWORDS ...], --exclude_keywords [EXCLUDE_KEYWORDS ...]
                        Keywords that, if present, will cause a line to be excluded.
  -l LOG_FILE, --log_file LOG_FILE
                        Output file to log rejected clusters (optional).
  --clean_words         Clean clusters based on a predefined dictionary of words.
  --all_species         If set, keeps the entire cluster if it contains at least one valid line.
  --selected_species    If set, removes clusters with only excluded lines but keeps clusters with mixed lines, removing only excluded lines.
```

**Explicaton supplémentaire:**
- -i
Ce paramètre permet de spécifier les rangs taxonomiques d'intérêt si l'on choisit de filtrer avec --all_species ou --selected_species.

- -e 
Ce paramètre permet d'exclure des mots qui ne seraient pas encore présents dans le dictionnaire clean_words_dict (voir ci-dessous).

- --clean_words
Le dictionnaire de mots utilisés pour le nettoyage des clusters est :
```
clean_words_dict = {
    'sp.': True,
    'Unassigned': True,
    'unknown': True,
}
```
Ces termes indiquent que les espèces ne sont pas encore bien annotées. On peut donc choisir de les supprimer afin d'obtenir uniquement des taxonomies propres.


**Méthodes de nettoyage des clusters**
- --all_species

Cette option supprime les clusters contenant uniquement des espèces non intéressantes.
Cependant, si un cluster contient à la fois des espèces d’intérêt et des espèces non pertinentes, il est conservé.

Cela permet d’évaluer la discrimination de nos espèces d’intérêt par rapport à celles capturées de manière non intentionnelle par les primers, sans fausser les statistiques avec des espèces non pertinentes. En effet, nous ne cherchons pas à discriminer les espèces qui ne nous intéressent pas, donc les clusters ne contenant que ces dernières ne nous concernent pas.

- --selected_species

Ce mode correspond à une situation où l'on travaille avec un échantillon contenant uniquement les espèces d’intérêt.
Dans ce cas, toutes les espèces non pertinentes sont supprimées.


##### Exemple avec nos data_test

Prenons l'exemple des Lactobacillus :

Tous les membres de ce groupe n’appartiennent plus au même genre qu’auparavant, où ils étaient tous classés sous le genre Lactobacillus. Avec la nouvelle classification, plusieurs genres ont été créés à partir de ce dernier.

Si je souhaite que mon fichier de statistiques ne représente que les Lactobacillus, je vais exécuter la commande suivante :

```
python ../script_treatment_ecopcr_result/taxonomic_discrimination_at_species.py -c cluster.txt -o stats.txt -i f__Lactobacillaceae g__Acetilactobacillus g__Agrilactobacillus g__Amylolactobacillus g__Apilactobacillus g__Bombilactobacillus g__Companilactobacillus g__Fructilactobacillus g__Furfurilactobacillus g__Lacticaseibacillus g__Lactiplantibacillus g__Lactobacillus g__Lapidilactobacillus g__Latilactobacillus g__Lentilactobacillus g__Levilactobacillus g__Ligilactobacillus g__Limosilactobacillus g__Liquorilactobacillus g__Loigolactobacillus g__Paucilactobacillus g__Philodulcilactobacillus g__Paralactobacillus g__Schleiferilactobacillus g__Secundilactobacillus g__Xylocopilactobacillus g__Holzapfelia g__Dellaglioa -l rejected_clusters.txt -r cluster_corrected.txt --clean_words --selected_species
```

Avec la liste de tous les genres que je souhaite conserver dans -i, ainsi que les paramètres --clean_words et --selected_species, je vais supprimer tout ce qui n'est pas Lactobacillus.

J'obtiendrai :

- Le fichier cluster_corrected.txt, qui ne contiendra que les espèces appartenant aux genres sélectionnés dans leur taxonomie.
- Le fichier rejected_clusters.txt, qui regroupera toutes les espèces ne faisant pas partie de ces genres.

Le fichier rejected_clusters.txt inclura donc toutes les espèces amplifiées par nos amorces qui ne sont pas des Lactobacillus.

NB: Il n'est pas obligatoire d'appliquer des inclusions et des exclusions spécifiques. Nous pouvons choisir de conserver tous les résultats sans modification.



##### Analyse des résultats du fichier stats.txt

Voici comment fonctionne le script :

Remarque : Nous travaillons sur la taxonomie complète. Si nous utilisons uniquement le taxid, les résultats peuvent varier, car deux taxid différents peuvent avoir la même taxonomie s'ils correspondent à des taxid de souche plutôt que d'espèce. Le script, cependant, utilise une taxonomie qui va jusqu'au niveau de l'espèce.

Le script évalue si une taxonomie est "bien discriminée" ou "mal discriminée" en examinant la diversité des taxonomies présentes dans chaque cluster.

Pour chaque cluster, le script extrait toutes les taxonomies associées aux séquences. Si toutes les séquences d'un cluster appartiennent à la même taxonomie (c'est-à-dire qu'il n'y a qu'une seule taxonomie unique dans le cluster), le cluster est considéré comme ayant une "bonne discrimination". La taxonomie est alors ajoutée à un ensemble de "taxonomies bien discriminées".

Si un cluster contient plusieurs taxonomies différentes (c'est-à-dire que plusieurs taxonomies sont présentes dans un même cluster), le cluster est considéré comme ayant une "mauvaise discrimination". Le script ajoute alors les taxonomies présentes dans ce cluster à un ensemble de "taxonomies mal discriminées".

Le script utilise des ensembles (set) pour stocker les taxonomies bien discriminées (good_discriminated_taxonomies) et mal discriminées (bad_discriminated_taxonomies).
Les ensembles sont des structures de données qui ne permettent pas de stocker des éléments en double. Si une taxonomie est ajoutée plusieurs fois à un ensemble, elle ne sera comptée qu'une seule fois.

Taxonomies ambiguës : Si une taxonomie apparaît à la fois dans des clusters bien discriminés et dans des clusters mal discriminés, elle est considérée comme "ambigüe" ou ayant des discriminations mixtes.


Résultat pour notre commande précédente:
```bash=
Number of unique taxonomies: 360 #Nombre total de taxonomies différentes dans les clusters
Number of well-discriminated taxonomies: 332
Percentage of well-discriminated taxonomies: 92.22%
Number of taxonomies both well and poorly discriminated: 47
Total number of filtered clusters: 633
Number of clusters with good discrimination: 582
Number of clusters with poor discrimination: 51
Percentage of clusters with good discrimination: 91.94%

# En dessous, nous avons les clusters mal discriminés ainsi que les espèces qui présentent à la fois une bonne et une mauvaise discrimination.
```

## 2. Comparaison avec la région V3V4 de l'ARNr 16S (Exemple)

Au cours de l'étape launch_format_ecopcr_result, un fichier de sortie nommé all_genome_name.txt est généré. Il contient tous les identifiants des espèces amplifiées lors de notre PCR in silico.

L'idée est d'évaluer si la région V3-V4 de l'ARNr 16S aurait été plus efficace pour discriminer ces espèces.

```
 cd 2_Amplicon_comparison/
```

### a. Extraction des séquences d'ARNr 16S dans la base de données
***Environnement conda TaxonMarker_main***

example de commande:
```
python amplicon_sequence_selected_extractor.py --genome_file ../1_EcoPCR/all_genome_name.txt --tsv_file $PATH_BD_16S/16S_with_taxonomy.tsv --fasta_file $PATH_BD_16S/16S.fna --output_file filtered_16S.fna
```

### b. Formatage du fasta obtenu avec obiconvert pour EcoPCR
***Environnement conda TaxonMarker_ecopcr***

example de commande
```
mkdir ecoPCR_db/
obiconvert --fasta filtered_16S.fna --ecopcrdb-output=ecoPCR_db/16S -t $PATH_NCBI_TAX_DUMB
```

### c. Lancement de la pcr in silico

example de commande
```
ecoPCR -d ecoPCR_db/16S CCTACGGGNGGCWGCAG GACTACHVGGGTATCTAATCC > 16S.ecopcr
```

### d. Analyse des résultats
***Environnement conda TaxonMarker_main***
Ce sont exactement les mêmes scripts que pour nos amorces:

example de commande
```
python ../script_treatment_ecopcr_result/format_ecopcr_result.py -o . -t name_seq_amplicon.txt 16S.ecopcr


python ../script_treatment_ecopcr_result/Launch_swarm.py -f all_modified.fna -s fichier_swarm.txt -o cluster.txt  -t 4 -a 1 -d 1 -vsearch output_vsearch_cluster.txt -taxo name_seq_amplicon.txt


python ../script_treatment_ecopcr_result/taxonomic_discrimination_at_species.py -c cluster.txt -o stats.txt -i f__Lactobacillaceae g__Acetilactobacillus g__Agrilactobacillus g__Amylolactobacillus g__Apilactobacillus g__Bombilactobacillus g__Companilactobacillus g__Fructilactobacillus g__Furfurilactobacillus g__Lacticaseibacillus g__Lactiplantibacillus g__Lactobacillus g__Lapidilactobacillus g__Latilactobacillus g__Lentilactobacillus g__Levilactobacillus g__Ligilactobacillus g__Limosilactobacillus g__Liquorilactobacillus g__Loigolactobacillus g__Paucilactobacillus g__Philodulcilactobacillus g__Paralactobacillus g__Schleiferilactobacillus g__Secundilactobacillus g__Xylocopilactobacillus g__Holzapfelia g__Dellaglioa -l rejected_clusters.txt -r cluster_corrected.txt --clean_words --selected_species
```

Résultats:

```
Number of unique taxonomies: 340
Number of well-discriminated taxonomies: 143
Percentage of well-discriminated taxonomies: 42.06%
Number of taxonomies both well and poorly discriminated: 26
Total number of filtered clusters: 272
Number of clusters with good discrimination: 201
Number of clusters with poor discrimination: 71
Percentage of clusters with good discrimination: 73.90%
```
