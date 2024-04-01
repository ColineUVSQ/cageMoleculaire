#include "utils_cage_moleculaire.h"
#include "graphe_cycles.h"
#include "lecture_molecule_sdf.h"
#include <dirent.h>

#define NB_MOL 81
#define NB_TAB 20

#define FILES_SMI "data/smi_files_reduit/"
#define FILES_DOT_GCYCLES "data/dot_files_reduit/graphes_cycles/"
#define FILES_DOT_GCOINS "data/dot_files_reduit/graphes_coins/"
#define FILES_PNG_GCYCLES "data/png_files_reduit/graphes_cycles/"
#define FILES_PNG_GCOINS "data/png_files_reduit/graphes_coins/"

#define RESULTS_ENS_CLIQUES "results/results_cliques_reduit.csv"
#define RESULTS_TYPE_CLIQUES "results/results_cliques_type_reduit.csv" 
#define RESULTS_DL_CLIQUES "results/results_clique_dl_reduit.csv"
#define RESULTS_LISTE_COINS "results/liste_coins_reduit.csv"
