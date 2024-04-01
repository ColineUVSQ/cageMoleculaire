#!/bin/bash

## creation des repertoires
echo "creation des repertoires"
for i in "data/smi_files_reduit/" "data/dot_files_reduit/" "data/dot_files_reduit/graphes_cycles/" "data/dot_files_reduit/graphes_coins" "data/png_files_reduit/" "data/png_files_reduit/graphes_cycles" "data/png_files_reduit/graphes_coins" "results/" 
do 
     if [ -d $i ]; then
    	echo "Le dossier existe ($i)"
     else
	echo "Le dossier n'existe pas ($i). Il va être créé"
	mkdir $i
     fi
done

#on cree les fichiers smi (nom_mol smiles) a partir des smiles du fichier d entree
echo "creation des fichiers smi et mol si pas encore fait"
python scripts/script_creer_petits_fic_smiles.py data/lot_cageV2.csv

## calcul des fichiers .mol (positions 3D et liaisons) a partir des fichiers .smi 
echo "fichiers mol: peut prendre un peu de temps"
bash scripts/script_smi_to_mol.sh
 
## on lance le programme
make run_cage ## calcul les graphes de cycles et graphes de coins du dataset lot_cageV2.csv
make clean

python scripts/script_compter_cliques.py ## compter les cliques de chaque type
bash scripts/genere_images_g_cycles.sh ## creer les fichiers png des graphes de cycles
bash scripts/genere_images_g_coins.sh ## creer les fichiers png des graphes de coins

rm scripts/genere_images_g_cycles.sh
rm scripts/genere_images_g_coins.sh


