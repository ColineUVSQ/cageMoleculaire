#ifndef NOM_D1
#define NOM_D1
#include "structure.h"


GRAPHE_CYCLE construction_graphe_cycles(struct molecule m);
void elimination_feuilles(struct molecule m);
graphemol conversion_mol_graphe(struct molecule m);
void probleme_memoire();
int *calcul_degre_mol( graphemol m);
int position_graphemol(graphemol g,int element);
graphemol modification_structure_mol(graphemol g, int* degre,int* deja_elimine);
int verif_mol(graphemol g, int* degre, int *deja_elimine);
graphemol suppression_aretes_mol(graphemol m, int i,int *degre);
void liberer_graphemol(graphemol g);
graphemol * ensemble_connexe_graphemol(graphemol m , int *deja_elimine);
struct liste_voisins*  construction_voisinage_graphemol(graphemol M);
int * sommets_atteignable_graphemol(int i,graphemol m,struct liste_voisins * v,int *deja_elimine);
void liberation_liste_voisins_graphemol(struct liste_voisins *v,graphemol m);
int nombre_isthmes_graphemol(graphemol g, int *deja_elimine);
int retirer_liaison_connexe(graphemol g, struct liaison l, int *deja_elimine);
graphemol enlever_une_arete(graphemol g, struct liaison l);
int est_connexe_graphemol(graphemol m, int *deja_elimine);
int existe_chaine_graphemol(int i, int j , graphemol m,struct liste_voisins * v);
graphemol ajouter_une_arete(graphemol g, struct liaison l);
isthmes *retrouver_tous_isthmes(graphemol *liste_connexe,isthmes *liste, int *deja_elimine);
graphemol enlever_tous_isthmes(graphemol g,int *deja_elimine,isthmes *liste_isthmes,int nb_isthmes);
int sommet_dans_basniveau(graphemol t,int p);
int *plus_court_chemin(int sommet1,int sommet2,graphemol m);
void liberer_memoire_voisins(struct liste_voisins *v,graphemol m);
int chemin_independant(int *chemin1, int *chemin2,struct liaison l);
cycles creer_un_cycle(graphemol m , int sommet, struct liaison l, int *chemin1,int *chemin2);
int verification_ajout_cycle(cycles *l, int nb_cycles , cycles c,graphemol m );
cycles *ajouter_un_cycle(cycles *liste, int nb_cycles , cycles c,graphemol m);
void liberer_un_cycle(cycles c);
int *concatener_deux_chemins(cycles c);
int position_de_arete(int sommet1, int sommet2,graphemol m);
int fonction_xor(int a , int b);
int *produit_xor_matrice(int *t, int *tab1, int *tab2,int nb);
int verification_egalite_tableaux(int *tab1,int *tab2 , int taille);
void ajouter_cycle_base(cycles c,graphemol m);
void arete_dans_cycle(struct molecule m);
void arete_dans_cycle_liste(struct molecule m);	
int sommets_commun ( cycles a, cycles b);
void nouvelle_arete(ARETE a);
int position_graphemol_arete(graphemol g, int sommet1,int sommet2);
int distance_inter_magma(cycles a , cycles b , struct molecule m);
int verification_LC(cycles a, cycles b,struct molecule m,int dist);
int nombre_de_cycles(int sommet1,int *chemin,int sommet2,graphemol g);
ARETE copier_arete(ARETE a , ARETE b);
void liberer_graphe_cycles( GRAPHE_CYCLE c);
int min( int a , int b);

void affichage_cycle(cycles c);
int compte_1_ligne(int** matrice, int nb_colonnes, int num_ligne);
int* copie_chemin(int* chemin_vieux);
int est_dans_liste_sommets(SOMMET_VR* liste_s, int nb_sommet_vr, int x);
int* init_t_zeros(int nb_liaisons);
int intersection_vide_chemins(int* chemin1, int* chemin2, int ext1, int ext2);
void libere_graphe_dr(struct graphe_dr* g_dr, int nb_atomes);
int** list_path(int x, int r, int* chemin, struct graphe_dr g_dr, int* nb_chemins_courants, int taille_chemin, int** chemins_courants, int l_A);
int nb_1_in_matrix(int** matrice, int i, int nb_col);
void obtenir_la_base_prototype(cycles *liste,int nb_cycles, graphemol m, struct graphe_dr* g_dr_i);
int *plus_court_chemin_precedant(int sommet1,int sommet2,graphemol m);
void trouver_prototypes_cycles_vismara(graphemol t);


#endif








