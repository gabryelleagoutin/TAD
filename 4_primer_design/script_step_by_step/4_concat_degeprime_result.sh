#!/bin/bash


# Définir le chemin du répertoire contenant les fichiers
directory="degeprime_result"

# Récupérer les préfixes des fichiers
prefixes=$(ls -1 "$directory" | grep -E ".*\.tsv$" | cut -d'_' -f1 | sort | uniq)

# Boucle sur chaque préfixe
for prefix in $prefixes; do
    # Concaténer tous les fichiers avec le même préfixe
    files=$(ls -1 "$directory"/"$prefix"* 2>/dev/null | grep -E ".*\.tsv$")
    if [ ! -z "$files" ]; then
        # Créer un nouveau fichier concaténé
        concatenated_file="$directory/concatenated_$prefix.tsv"
        # Concaténer tous les fichiers en enlevant l'en-tête sauf pour le premier
        cat $(echo "$files" | head -n1) > "$concatenated_file"
        for file in $(echo "$files" | tail -n+2); do
            tail -n +2 "$file" >> "$concatenated_file"
        done
    fi
done
