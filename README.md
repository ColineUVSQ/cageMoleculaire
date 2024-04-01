## Construction graphe de cycles et graphe de coins pour la recherche de cages moléculaires

Le programme prend en entrée un fichier .csv contenant des molécules et leur SMILES (fichier data/lot_cageV2.csv par défaut, 81 molécules). A partir de ce fichier, le programme effectue les étapes suivantes:

- Génération des structures 3D des molécules (coordonnées des atomes et liaisons entre atomes) avec le package OpenBabel et stockage dans des fichiers .mol (pas forcément nécessaire, il faudrait avoir un programme qui passe directement du SMILES au graphe moléculaire sinon)
- Calcul du graphe moléculaire de chaque molécule à partir du fichier .mol (Cette étape prend quelques minutes)
- Calcul du graphe de cycles de chaque molécule, comme défini dans [1], avec comme ensemble de cycles l'union des bases de cycles [2]. 
- Calcul du graphe de coins de chaque molécule, et classification des coins selon leurs types (voir fichier data/Détection_de_cages_moléculaires_à_l_aide_de_graphes_de_cycles.pdf).


NB: Dans le fichier lot_cageV2.csv, les SMILES utilisés sont ceux de la colonne "smiles ring system", correspondant à la zone de la molécule pouvant ou non former une cage.

### Pour utiliser le programme
Exécuter le fichier **script.sh** permet de réaliser l'ensemble des étapes.
En ligne de commande:
```
    chmod +x script.sh
    ./script.sh
```

### Sortie du programme
Le programme génère en sortie plusieurs fichiers csv (dossier results):
D'abord des fichiers intermédiaires:
- *results_cliques_reduit.csv*: stocke l'ensemble des cliques de taille 3 du graphe de cycle trouvées par molécule
- *results_cliques_type_reduit.csv*: stocke l'ensemble des types de clique de taille 3 du graphe de cycle trouvées par molécule
- *results_clique_dl_reduit.csv*: stocke les cliques de taille 3 du graphe de cycle comprenant des double liaisons (à la jonction de 2 cycles ou dans un cycle seulement)

Et des fichiers directement interprétables:
- *results_cliques_type_compte_reduit.csv*: stocke le nombre de cliques de chaque type par molécule
- *liste_coins_reduit.csv*: stocke les poids des sommets et les poids des arêtes de chaque coin par molécule

Des fichiers .png pour visualiser les graphes de cycles et les graphes de coins sont également générés (dossier data).

### Références
[1] Stefi Nouleho Ilemo, Dominique Barth, Olivier David, Franck Quessette, Marc-Antoine Weisser,
and Dimitri Watel. Improving graphs of cycles approach to structural similarity of molecules.
PLOS ONE, 14(12):e0226680, December 2019. Publisher: Public Library of Science.

[2] Philippe Vismara. Union of all the Minimum Cycle Bases of a Graph. The Electronic Journal
of Combinatorics, pages R9–R9, January 1997.
