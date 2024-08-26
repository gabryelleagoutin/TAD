# TaxonMarker

[TOC]

## Présentation
TaxonMarker est un outil développé en Python et Bash, conçu pour identifier des marqueurs génétiques spécifiques à un rang taxonomique donné. Lors de l'analyse de séquençage métabarcoding, le gène de l'ARNr 16S est couramment utilisé, bien qu'il présente des limites, notamment sa présence en copies multiples et sa faible résolution au niveau des organismes étroitement apparentés. Le but de TaxonMarker est de selectionner des gènes cibles potentiellement plus efficaces que le 16S pour discriminer les espèces. En entrée, il suffit de fournir l'identifiant taxonomique souhaité (par exemple, Actinomycetes, taxid: 1760).

## Installation 

```bash=
git clone https://github.com/GTG1988A/TaxonMarker.git
```
TODO : Environnement conda ? 

## 1. formatting file bacterian OG

Cette étape doit être effectuée uniquement par le programmeur, à chaque mise à jour microbienne d'OrthoDB. Elle permet de récupérer uniquement les groupes orthologues bactériens dans la base de données OrthoDB.

```bash=
cd TaxonMarker/1_formatting_file_bacterian_OG
python formatting_bacterial_orthologue_file_final.py -o ../Orthodb/odb11v0_OG2genes.tab -b only_line_bacteria.txt -u uniq_og_ids.txt -s ../Orthodb/odb11v0_level2species.tab -g ../Orthodb/odb11v0_OGs.tab -f Bacterial_OG.tab
```

## 2. search taxid and monocopy calculation

Récupération des OGs contenant le rang taxonomique sélectionné. On extrait les identifiants des espèces issues de ce rang taxonomique, puis on conserve les OGs qui contiennent au moins une de ces espèces. Cette liste d'espèces est appelée identifiers_with_searchID_in_taxonomy.

Cette liste identifiers_with_searchID_in_taxonomy est également utilisée pour vérifier combien de taxid dans la colonne "taxid" appartiennent à cette liste. Cela nous permet de déterminer le nombre d'espèces appartenant à notre rang souhaité dans l'OG, ainsi que le pourcentage par rapport au nombre total d'espèces.

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
python search_taxid_and_monocopy_and_percentage_calculation.py -i ../1_formatting_file_bacterian_OG/Bacterial_OG.tab -f /BD_TaxonMarker/Orthodb/odb11v0_species.tab -l /BD_TaxonMarker/Orthodb/odb11v0_level2species.tab -s 1578 -o OG_1578.tab
```
### fonctionnement de OG_selection.sh et Astuces
Le script peut trier selon trois critères :

> -p <percent_single_copy> (obligatoire, mettre à 0 si aucune sélection n'est souhaitée) 
> [-c <TargetSpecies_Count>] Nombre minimum d'espèces (car un OG contenant seulement deux espèces à 100 % de notre rang est moins intéressant qu'un OG en contenant 100) 
> [-t <TargetSpecies_Percentage>] Pourcentage minimum pour n'inclure que les OGs correspondant potentiellement à notre rang d'intérêt

On peut déjà savoir combien d'espèces descendent de notre rang taxonomique en utilisant la commande suivante (elle est indiquée dans le --help du script).
```bash=
$ grep -w "1578" /BD_TaxonMarker/Orthodb/odb11v0_levels.tab
1578    Lactobacillus   542691  14672   265
```

On peut donc essayer d'aller à fond sur les paramètres du genre:

```bash=
./OG_selection.sh -p  100 -c 265 -t 100 -o OG_1578_selected.tab OG_1578.tab
```

puis réduire progressivement jusqu'à obtenir un résultat satisfaisant.


## 3. fasta recovery

Téléchargement des séquences des gènes au format acide nucléique. Pour cela, nous utilisons deux API. Ce script récupère les OGs sélectionnés à l'étape 2. Il extrait l'ID protéique de chaque protéine dans l'OG, puis utilise l'API d'OrthoDB pour obtenir l'ID EMBL du CDS. Cet identifiant est ensuite utilisé pour télécharger la séquence nucléique au format FASTA via l'API EMBL. Si aucun ID n'est trouvé, cela sera indiqué dans les logs.

*Exemple de commande:*
```bash!
python fastas_recovery.py ../2_search_taxid_and_monocopy_calculation/OG_1578_selected.tab
```

Le script ajoute aussi dans le tableau le nombre de séquence contenues dans le fasta de l'OG.

## 4. primer design

Cette étapes en contient plusieurs:

### a. Installation de Degeprime 

```bash
git clone https://github.com/EnvGen/DEGEPRIME.git
```

### b. Alignement
Lancement de l'alignement des séquences de tout les OGs avec Clustal Omega

*Exemple de commande:*
```bash!
mkdir -p alignment/

for files in ../3_fasta_recovery/*.fa; do
    file_name=$(basename "$fichier")
    clustalo -i "$files" -o "alignment/$file_name" --threads=4
done
```
### c. Lancement de trim de Degeprime 

*NB: Voir + à propos de trim dans le github de degeprime.*
*Exemple de commande:*
```bash!
for files in alignment/*.fa; do
    file_name=$(basename "$fichier" .fa)
    perl DEGEPRIME/TrimAlignment.pl -i "$files" -min 0.9 -o "alignment/trimmed_$file_name.fna"
done

```
### d. Génération de toutes les amorces possibles

Toutes les amorces possibles sont générées en spécifiant les paramètres de longueur minimale et maximale du primer, ainsi que les valeurs de  dégénérescences.

La commande ci-dessous permet de generer un fichier sarray afin de lancer toutes les combinaisons de paramètres possible.


*Exemple de commande:*


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


NB: La dégénérescence correspond au nombre total de combinaisons possibles de nucléotides qu'un primer dégénéré peut former. Par exemple, pour le primer ATCS, où S représente G ou C, les combinaisons possibles sont ATCG et ATCC. La dégénérescence est donc égale à 2. Les valeurs possibles de dégénérescence peuvent être 2, 4, 8, 16, 32, 64, 96, etc. Il est possible d'atteindre des valeurs plus élevées, mais une dégénérescence de 96 est déjà extrêmement peu stringeante.

### e. Concaténation des résultats

On concatène les résultats pour avoir un fichier par OGs

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
### f. Ajout de métrique sur les primers

Nous lançons un programme qui ajoutera des métriques sur les primers 

Les statistiques ajoutées sont: température minimale et maximale, pourcentage de GC, présence d'un GC clamp, pourcentage de number matching, si le primer fini par un T, si il est complémentaire à lui même.

Nous spécifierons les températures maximales et minimales souhaitées, ainsi que le pourcentage minimum de correspondance (= nombre de séquences où une amorce est trouvée de manière identique). Les primers ne répondant pas à ces critères seront supprimés.

Nous obtiendrons un fichier TSV contenant les primers correspondant à nos critères, avec les nouvelles statistiques ajoutées dans le tableau.

*Exemple de commande:*
```bash!
python process_primers_stat.py -i degeprime_result/concatenated_* -og ../3_fasta_recovery/updated_OG_selected_1578.tab -o result_stat_primers -nm 80 -tm_max 65 -tm_min 54
```

### g. Création du tableau des couples de primers.

Nous lançons un programme qui générera les couples de primers possibles. Nous spécifierons la taille minimale et maximale de l'amplicon entre les deux primers. Les couples ne respectant pas ces critères seront supprimés. De plus, les couples avec un écart de température supérieur à 5 degrés entre le primer forward et le primer reverse seront également éliminés.

*Exemple de commande:*
```bash!
#!/bin/bash
for file in result_stat_primers/concatenated_*; do echo "module load devel/python/Python-3.11.1;python couple_primer.py -i $file -f alignment/ --amplicon_min_size 150 --amplicon_max_size 590" >> tab_couple.sarray;done

sarray --meme=200G tab_couple.sarray
```

Chaque couple de primers se verra attribuer un score total, calculé en additionnant le score pour le number matching et celui pour la taille de l'amplicon. Nous avons jugé que ces paramètres étaient les plus importants après avoir effectué la sélection basée sur les températures et les tailles d'amplicon.


### h. Concaténation des résultats et séléction des meilleurs couples.

On concatène les résultats avec concatenate_sort_result.sh et on selectionne les couples qui nous intéressent.

Nous sélectionnerons les couples ayant les trois couples meilleurs scores. Par la suite, nous pourrons essayer de choisir ceux ayant la dégénérescence la plus faible, un pourcentage de GC de 50 % ou un GC clamp à la fin, etc. Cependant, il est important de noter que le primer parfait n'existe pas nécessairement.

Un script nommé script_js_html_two_tab.py prend en entrée la liste des primers sélectionnés. Il génère un histogramme montrant la répartition de la dégénérescence. En effet, il est préférable que la dégénérescence soit plus concentrée  au début du primer, plutôt qu'à la fin. Cet histogramme permet donc une meilleure visualisation de cette répartition.

TODO: graphe sur la teneur en GC? Position?




NB : Si vous exécutez ce programme sur un cluster de calcul, vous pouvez utiliser le script primer_pipeline.sh en ajustant les paramètres à l'intérieur. Cela permet de lancer toutes les étapes en une seule fois.

## 5. EcoPCR

### a. Lancement de EcoPCR
Lancement de EcoPCR sur une base de données formatée pour l'outil (voir la section "Base de données pour TaxonMarker").

Cela nous permet de tester nos amorces in silico. Nous pouvons vérifier si elles capturent correctement les cibles souhaitées et si elles capturent également d'autres séquences non prévues.

Les EcoPCR sont généralement lancées une par une avec les primers sélectionnés. (À discuter si nous mettrons en place un processus automatisé pour lancer plusieurs EcoPCR simultanément.)

Pour lancer une EcoPCR, utilisez le script 1_commande_ecopcr.sh.
```bash=
./1_ecopcr_command.sh
```

Il ressemble a ça:

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

Vous pouvez spécifier un nom pour le fichier de sortie en définissant s_array_file. Une fois la première exécution réalisée, il est possible de lancer d'autres analyses.

Il est également nécessaire de fournir les informations suivantes :

    Le chemin de la banque EcoPCR (fourni par nous)
    Le nom du gène cible ou de son OG
    Les amorces forward et reverse sélectionnées

Le script génère alors la commande EcoPCR que vous pouvez exécuter comme suit :

```bash=
sarray --mem=200G ecopcr_commands.sarray
```

Plusieurs scripts sont a disposition pour traiter l'information:

### b. Analyse des résultats d'EcoPCR

#### 1_launch_format_ecopcr_result.sh

Ce script permet de convertir les résultats d'EcoPCR au format FASTA, d'ajouter les informations de taxonomie, puis de lancer l'outil vsearch pour la déréplication des séquences. Il supprime également toutes les séquences contenant des acides aminés, ce qui est nécessaire pour que SWARM fonctionne correctement.

*Exemple de commande:*

```bash=
python /path/Lactobacillus/TaxonMarker/script_treatment_ecopcr_result/format_ecopcr_result.py -o . -t /path/BD_TaxonMarker/BD_ecoPCR/name_seq_with_taxo.txt primers/*ecopcr

```

Il prend en argument le fichier name_seq_with_taxo.txt présent dans la base de données EcoPCR_DB pour TaxonMarker (voir la section "Base de données"). Ce fichier contient les noms des séquences au format FASTA accompagnées des informations de taxonomie.

Le fichier de sortie final est nommé all_modified.fna.


#### 2_launch_swarm.sh

Nous utilisons le clustering avec l'outil Swarm pour mesurer la discrimination des espèces. Swarm est un outil hautement discriminant, robuste et rapide, qui regroupe les amplicons ayant des ressemblances élevées entre eux. Nous partons du principe que deux amplicons regroupés dans le même cluster proviennent de la même espèce.

*Exemple de commande:*
```bash!
python /PATH/TaxonMarker/script_treatment_ecopcr_result/Launch_swarm.py -f all_modified.fna -s fichier_swarm.txt -o cluster.txt  -t 4 -a 1 -d 1
```

#### 3_launch_stat_swarm.sh

Ce script nous permet d'obtenir des mesures de discrimination des espèces spécifiquement pour notre rang taxonomique. Il offre la possibilité d'inclure ou d'exclure certains termes dans la taxonomie, afin de générer un fichier statistique propre, limité uniquement à notre rang d'intérêt.

##### Astuce au travers d'un cas concret:
Prenons l'exemple des Lactobacillus:

Tous les membres de ce groupe n’appartiennent plus au même genre qu’auparavant, où ils étaient tous classés sous le genre Lactobacillus. Depuis la nouvelle classification, plusieurs genres ont été créés à partir de ce genre. Cependant, ils partagent toujours la même famille, Lactobacillaceae.


Si je lance mon script en incluant uniquement la famille Lactobacillaceae donc en mettant
-i f__Lactobacillaceae

J'aurai le fichier cluster_corrected.txt qui contient uniquement des espèces appartenant à la famille Lactobacillaceae dans leur taxonomie, et le fichier rejected_clusters.txt avec celles qui n'appartiennent pas à cette famille.

Cependant, en examinant les résultats, nous constatons que la famille Lactobacillaceae comprend non seulement des Lactobacillus, mais aussi d'autres genres que nous souhaitons exclure. Nous dresserons alors une liste de tous les genres ou espèces à retirer:

g__Leuconostoc g__Pediococcus g__Oenococcus g__Unassigned sp. g__Weissella g__Dellaglioa g__Periweissella g__Convivina


On relance donc notre script avec l'exclusion de ces mots là

```bash=
python stats_report_taxo2.py -c cluster.txt -t all_modified.fna -o stats.txt -i f__Lactobacillaceae -e g__Leuconostoc g__Pediococcus g__Oenococcus g__Unassigned sp. g__Weissella g__Dellaglioa g__Periweissella g__Convivina -l rejected_clusters.txt -r cluster_corrected.txt

```

Nous obtiendrons alors le fichier cluster_corrected.txt contenant uniquement des espèces du genre Lactobacillus. Les statistiques seront calculées à partir de ce fichier.

Le fichier rejected_clusters.txt comprendra toutes les espèces amplifiées par nos amorces qui ne sont pas des Lactobacillus.

Remarque : Nous ne sommes pas obligés d'appliquer des inclusions et des exclusions spécifiques. Nous pouvons choisir de conserver tous les résultats sans modification.

NB : Par exemple pour la classe "Actinomycete" il n'y aura pas besoin de mots d'exclusion. 


##### Analyse des résultats du fichier stats.txt

Voici comment fonctionne le script :

Remarque : Nous travaillons sur la taxonomie complète. Si nous utilisons uniquement le taxid, les résultats peuvent varier, car deux taxid différents peuvent avoir la même taxonomie s'ils correspondent à des taxid de souche plutôt que d'espèce. Le script, cependant, utilise une taxonomie qui va jusqu'au niveau de l'espèce.

Le script évalue si une taxonomie est "bien discriminée" ou "mal discriminée" en examinant la diversité des taxonomies présentes dans chaque cluster.

Pour chaque cluster, le script extrait toutes les taxonomies associées aux séquences. Si toutes les séquences d'un cluster appartiennent à la même taxonomie (c'est-à-dire qu'il n'y a qu'une seule taxonomie unique dans le cluster), le cluster est considéré comme ayant une "bonne discrimination". La taxonomie est alors ajoutée à un ensemble de "taxonomies bien discriminées".

Si un cluster contient plusieurs taxonomies différentes (c'est-à-dire que plusieurs taxonomies sont présentes dans un même cluster), le cluster est considéré comme ayant une "mauvaise discrimination". Le script ajoute alors les taxonomies présentes dans ce cluster à un ensemble de "taxonomies mal discriminées".

Le script utilise des ensembles (set) pour stocker les taxonomies bien discriminées (good_discriminated_taxonomies) et mal discriminées (bad_discriminated_taxonomies).
Les ensembles sont des structures de données qui ne permettent pas de stocker des éléments en double. Si une taxonomie est ajoutée plusieurs fois à un ensemble, elle ne sera comptée qu'une seule fois.

Taxonomies ambiguës : Si une taxonomie apparaît à la fois dans des clusters bien discriminés et dans des clusters mal discriminés, elle est considérée comme "ambigüe" ou ayant des discriminations mixtes.


```bash=
Number of unique taxonomies: 291 #Nombre total de taxonomies différentes dans les clusters
Number of good discriminated taxonomies: 263
Percentage of good discriminated taxonomies: 90.38% 
Number of taxonomies both good and bad discriminated: 9
Total number of filtered clusters: 324 
Number of clusters with good discrimination: 304
Number of clusters with bad discrimination: 20
Percentage of clusters with good discrimination: 93.83%

# En dessous, nous avons les clusters mal discriminés ainsi que les espèces qui présentent à la fois une bonne et une mauvaise discrimination.
```

## 6. genome selection

Nous avons donc les informations sur la discrimination de nos gènes cibles pour le rang taxonomique sélectionné. Mais est-ce que le gène 16S n’est pas finalement meilleur ? Pour le savoir, nous devons effectuer une comparaison.

Lors de l'EcoPCR, nous testons nos amorces sur une base de données microbienne complète (voir la section sur la base de données). La première étape consiste à filtrer cette base de données pour conserver uniquement les espèces correspondant à notre rang taxonomique d'intérêt. Nous sélectionnons les génomes dont la taxonomie contient les critères pertinents (par exemple, en utilisant les mots d'inclusion et d'exclusion).

Cela nous permettra de rechercher, parmi ces génomes, ceux qui possèdent une séquence 16S. Ainsi, nous pourrons comparer correctement notre analyse du gène cible avec celle du 16S, puisque nous aurons sélectionné les génomes d'intérêt dans la même base de données.

*Exemple de commande:*
Par exemple, pour les Lactobacillus, nous lancerions la recherche sur notre base de données microbienne en utilisant les critères suivants :

```bash=
python filter_seq_with_keywords.py -s /home/gagoutin/work/BD_TaxonMarker/BD_ecoPCR/name_seq_with_taxo.txt -o filtered_seq_with_taxo.txt -i f__Lactobacillaceae -e g__Leuconostoc g__Pediococcus g__Oenococcus g__Unassigned sp. g__Weissella g__Dellaglioa g__Periweissella g__Convivina
```


TODO: Ce script nous sortira dans le log le nombre d'espèces total dans notre base de données. Cette information est intéressante, car on peut la comparer avec "Number of unique taxonomies:" de notre fichier stats.txt obtenu juste avant. Cela permet de voir combien on en attrape par rapport au total. 

NB pour les lacto:
```bash=
gagoutin@genobioinfo1 ~/work/Lactobacillus/TaxonMarker/6_genome_selection $ cut -d" " -f3 filtered_seq_with_taxo.txt | sort -u | wc -l
365
```

Le script produira également le fichier genome_name_selected.txt, qui contient la liste des noms des génomes sélectionnés. Ce fichier nous permettra de filtrer les génomes d'intérêt dans notre banque de séquences 16S, afin de conserver uniquement les séquences 16S correspondant aux rang souhaités (ex: Lactobacillus.)


## 7. 16S comparison
Nous commençons par lancer le script 16S_sequence_selected_extractor.py, qui prend en entrée notre liste de génomes précédemment obtenue. Ce script extrait uniquement les séquences 16S correspondant au rang taxonomique sélectionné.

```bash=
python 16S_sequence_selected_extractor.py --genome_file ../6_genome_selection/genome_name_selected.txt --tsv_file /PATH/BD_TaxonMarker/BD_16S/16S_with_taxonomy.tsv --fasta_file /PATH/BD_TaxonMarker/BD_16S/16S.fna --output_file filtered_16S.fna
```

Nous aurons plusieurs sorties :

    Le fichier FASTA filtré : filtered_16S.fna
    name_seq_16S.txt : les informations sur la taxonomie de chaque séquence
    process.log : où nous trouverons à la fin le "Number of unique taxonomies", qui indique le nombre d'espèces possédant un 16S

Par exemple, si nous avons 365 taxonomies différentes pour les Lactobacillus, nous obtenons 353 espèces avec un 16S.

Ensuite, nous formatons le fichier FASTA obtenu, en conservant uniquement les séquences 16S du rang sélectionné, à l'aide de la commande obiconvert.

```bash=
obiconvert --fasta filtered_16S.fna --ecopcrdb-output=ecoPCR_db/16S -t /PATH/BD_TaxonMarker/ncbi_tax_dumb/2024-15-04

```

Puis on lance ecopcr
NB : Mettre les amorces 16S couramment utilisé dans le laboratoire. 
```bash=
ecoPCR -d ecoPCR_db/16S CCTACGGGNGGCWGCAG GACTACHVGGGTATCTAATCC > 16S.ecopcr
```


### Traitement des résultats
On traite en suite les résultats comme précedemment:

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

NB : pas besoin de mettre des inclusions ou exclusion car on a filtrer avant.


#### Exemple analyse résultat pour Lactobacillus

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

Les clusters semblent plutôt bien discriminés, mais nous constatons que de nombreuses espèces identiques sont réparties sur plusieurs clusters. En réalité, les taxonomies bien discriminées ne correspondent qu'à 61 %.


### Comparaison avec species_comparison_venn.py

```bash!
python species_comparison_venn.py --taxo_target_gene ../5_EcoPCR/uniq_taxo.txt --taxo_16S_gene uniq_taxo.txt --good_discrimination_target_gene ../5_EcoPCR/uniq_taxo_good_discriminated.txt --good_discrimination_16S uniq_taxo_good_discriminated.txt --output venn_diagram.png
```


# Les différentes Database pour TaxonMarker
## OrthoDB

Base de données maintenue par l’Université de Genève et le SIB (Swiss institue of Bioinformatics)

Elle contient des groupes de gènes orthologues de tout le règne vivant
 = Tous les descendants d'un gène unique particulier du dernier ancêtre commun de ces espèces

C'est la base de données qui sert à l'outil TaxonMarker pour recuperer des gènes cibles. 

liens pour télécharger OrthoDB:

```bash
wget -r https://data.orthodb.org/v11/download/
```


## 16S Database 
### Description de la base de données utilisées (identique pour EcoPCR)
Nous allons extraire les 16S a partir d'une base de données contenant les Génomes microbiens NCBI (fasta et gff) contenus dans les bases de données RefSeq et GenBank avec la taxonomie GTDB avec l'outil genome_updater.

Cette bases provient de ce téléchargement 
```bash
genome_updater.sh -g "archaea,bacteria" -d "refseq,genbank" -M "gtdb" -f "genomic.gff.gz,genomic.fna.gz" -a -o downloads -N -t 2 -V -Z -R 6 -L curl
```
avec cette outils: 
https://github.com/pirovc/genome_updater
Bash script to download/update snapshots of files from NCBI genomes repository (refseq/genbank) with track of changes and without redundancy

Elle est conservé au cluster Genotoul de Toulouse, toute les informations sont disponible ici:

```bash 
/work/bank2/users_banks/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/README
```

### Extraction des 16S

Nous avons besoin des scripts
- regions_analysis_fct.py
- identify_rrna_gene.py
- extract_16S23S.py

le scripts extract_16S23S.py a besoin de deux fichiers
python extract_16S23S.py -a genome_dirs.txt -t genome_files_with_taxid.txt

genomes_dirs.txt contient la localisation de chaque fichier fna et gff des génomes
le fichier genome_files_with_taxid.txt contient le chemin jusqu'au nom du fichiers des génomes et en plus le nom du génomes et son taxid 

```bash=
# $ head genome_files_with_taxid.txt
/home/gagoutin/work/base_données//NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/GCA/018/630/415/GCA_018630415.1_ASM1863041v1_genomic.fna.gz GCA_018630415.1 1280
/home/gagoutin/work/base_données//NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/GCA/018/630/395/GCA_018630395.1_ASM1863039v1_genomic.fna.gz GCA_018630395.1 1280

# $ head genome_dirs.txt
/home/gagoutin/work/base_données/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/GCA/018/630/415
/home/gagoutin/work/base_données/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/2024-05-15_11-08-32/files/GCA/018/630/395

```

pour avoir ces fichiers, il faut d'abord lancer les scripts:
./1_find_dirs.sh
et 2_find_name.py en prenant bien soin d'ajuster les chemins de la base de données utilisées.

Une fois le script d'extraction lancé, nous obtenons deux fichiers résultants :

    16S.tsv : contenant les informations sur chaque région 16S extraite.
    16S.fna : contenant les séquences des 16S extraits.

Je vais également ajouter un fichier 16S_with_taxonomy.tsv qui ajoute une colonne de taxonomie au fichier 16S.tsv. Ce fichier nous permettra de sélectionner les taxonomies qui nous intéressent lors de la comparaison entre le 16S et un gène cible. Par exemple, pour les Lactobacillus, nous pourrons filtrer pour ne conserver que ces espèces et exclure les autres.

Pour obtenir ces taxonomies, j'utilise TaxonKit.

```bash=
cut -f 12 16S.tsv  | tail -n +2 | taxonkit lineage --data-dir ../ncbi_tax_dumb/2024-15-04/ | taxonkit reformat --data-dir ../ncbi_tax_dumb/2024-15-04/  -f "{k};{p};{c};{o};{f};{g};{s}" -r "Unassigned" -P | cut -f 1,3-10  > 16S.taxo

```

Lorsque la colonne est vide car pas de taxonomie j'ajoute un Unassigned:
```bash=
awk 'BEGIN {FS=OFS="\t"} {if ($2 == "") $2 = "k__Unassigned;p__Unassigned;c__Unassigned;o__Unassigned;f__Unassigned;g__Unassigned;s__Unassigned"; print}' 16S.taxo > 16S_taxo_modified.txt

```

J'enlève les espaces indésirable: 

```bash=
 sed -i 's/ /_/g' 16S_taxo_modified.txt
```

Je conserve uniquement la deuxieme colonne et j'y ajoute l'header "taxonomy"

```bash=
cut -f2 16S_taxo_modified.txt | sed  '1s/^/taxonomy\n/' > col2.txt
```

Je m'assure que les fichiers soient bien encodés correctement:

```bash=
dos2unix col2.txt 16S.tsv
```
/!\ Il est important d'effectuer cette vérification, car si le fichier TSV est encodé en "ASCII text, with CRLF line terminators" au lieu de simplement "ASCII text", la commande paste ne fonctionnera pas correctement. Vous pouvez vérifier l'encodage du fichier en utilisant la commande suivante :
```bash!
file 16S.tsv
```

et maintenant je colle la colonne de la taxonomie dans le fichier tsv:
```bash=
paste 16S.tsv col2.txt > 16S_with_taxonomy.tsv
```
## ncbi tax_dump

Cela nous permets d'avoir la taxonomie du ncbi.
toute les infos se trouvent ici:
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/

On le télécharge comme ça:
```bash=
wget -r https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/

```
## EcoPCR DataBase
### Description de la base de données utilisées (identique pour 16S)
Nous formatons une base de données contenant les Génomes microbiens NCBI (fasta et gff) contenus dans les bases de données RefSeq et GenBank avec la taxonomie GTDB avec l'outil genome_updater.

Cette bases provient de ce téléchargement 
```bash
genome_updater.sh -g "archaea,bacteria" -d "refseq,genbank" -M "gtdb" -f "genomic.gff.gz,genomic.fna.gz" -a -o downloads -N -t 2 -V -Z -R 6 -L curl
```
avec cette outils: 
https://github.com/pirovc/genome_updater
Bash script to download/update snapshots of files from NCBI genomes repository (refseq/genbank) with track of changes and without redundancy

Elle est conservé au cluster Genotoul de Toulouse, toute les informations sont disponible ici:
```bash 
/work/bank2/users_banks/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/README
```

Elle est formater comme ça:
### Formatter la base de données pour EcoPCR

#### 1. Récupérer l'assembly summary
```bash
cp  /work/bank2/users_banks/NCBI_refseq_genbank_archea_bacteria/NCBI_refseq_genbank_archea_bacteria_05_2024/downloads/assembly_summary.txt .
```

#### 2. Compléter avec l'header
```bash
cat header.txt
assembly_accession      bioproject      biosample       wgs_master      refseq_category taxid   species_taxid   organism_name   infraspecific_name      isolate version_status  assembly_level  release_type genome_rep       seq_rel_date    asm_name        asm_submitter   gbrs_paired_asm paired_asm_comp ftp_path        excluded_from_refseq    relation_to_type_material       asm_not_live_date       assembly_typegroup    genome_size     genome_size_ungapped    gc_percent      replicon_count  scaffold_count  annotation_provider     annotation_name annotation_date total_gene_count        protein_coding_gene_count    non_coding_gene_count    pubmed_id

cat header.txt assembly_summary.txt> assembly_summary_with_header.txt
```

#### 3. Séparation en morceaux pour aller plus vite en raison des très nombreux génomes.

```bash
nb_chunk=500
mkdir assembly_chunks
cut -f1 assembly_summary_with_header.txt > ref_and_rep_assemblies_acc_list.txt
split -a 3 -d --number=l/$nb_chunk ref_and_rep_assemblies_acc_list.txt assembly_chunks/ref_and_rep_assemblies.acc_list_chunk.
```


#### 4. La commande sarray pour formater la base de données
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
#### 5. Lancement du sarray

```bash
sarray --mem=200G build_assemblies_ecopcr_db.sarray
```
#### 6. Création de name_seq avec et sans taxonomie. 
Cela me permet d'avoir le contenu de l'espèce dans tous mes génomes et de conserver cette information pour des comparaisons futures.

Récupération du nom du génome à l'aide du taxid. Cela correspond au nom des séquences dans les fastas :
```bash
zcat *fna.gz | grep ">" | cut -d' ' -f1-2 > name_seq_without_taxo.txt
```

Utilisation du kit Taxon, qui prend le taxid en entrée et produit la taxonomie jusqu'au niveau de l'espèce :
```bash
cut -d'=' -f2 name_seq_without_taxo.txt | cut -d';' -f1 | taxonkit lineage –data-dir ../rpob/complete_and_chromosome/ncbi_tax_dumb/2024-15-04/ | taxonkit reformat –data-dir ../rpob/complete_and_chromosome/ncbi_tax_dumb/2024-15-04/ -f "{k};{p};{c};{o};{f};{g};{s}" -r "Unassigned" -P | cut -f 1,3-10 > name_seq.taxo
```

Supprimez tous ceux dont l'ID est introuvable, sinon obiconvert ne les prendra pas. Vous devez également supprimer les espaces dans les noms :

```bash
awk '{if (NF > 1 && $2 != "") print $0}' name_seq.taxo | sed 's/ /_/g' > filtered_name_seq.taxo
```

Lancer le script pour ajouter le taxo au nom de la séquence :
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



