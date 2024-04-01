#include "graphe_cycles.h"

struct mol_carac{
	char* name; // nom de la molécule
	int classe; // 1:cage 2:precage ou 3:non cage
};

typedef struct mol_carac MOL_CARAC;

MOL_CARAC *classification;

struct arete
{
	int id1;
	int id2;
	int poids;
};
typedef struct arete ARETE_COIN;

struct sommet
{
	int id;
	int* ids; // ids des sommets correspondants dans la clique du graphe de cycles
	int* poids; // poids des 3 cycles de la clique
	int* liaisons_communs; // nb de liaisons en commun entre les 3 cycles

};
typedef struct sommet SOMMET_COIN;

struct graphe_coin
{
	int nb_sommets;
	int nb_aretes;
	ARETE_COIN *liste_aretes;
	SOMMET_COIN *liste_sommets;
};

typedef struct graphe_coin GRAPHE_COIN;

// pile pour convertir smiles en molecular graph
typedef struct stack
{
   struct stack *prev;
   struct stack *next;
   int num;
} stack_s;


stack_s* init_stack(); // initialiser la pile
void stack_push(stack_s ** pp_stack, int num); // empiler
int stack_pop(stack_s ** pp_stack); // dépiler
void stack_delete(stack_s ** pp_stack); // supprimer la pile


void init_cage_non_cage();

// Returns a 2 dimensions tabular containing nb_lignes and nb_col with -1 everywhere
int** init_tab_deux_dimensions(int nb_lignes, int nb_col);

// Return the vertex a in which are copied the values of fields of the vertex b
SOMMET copier_sommet(SOMMET a, SOMMET b);

// Returns a tab with all the degrees of the vertices in the graph of cycles cy
int *calcul_degre(GRAPHE_CYCLE cy);

// Function that returns 1 if (x,y) is an edge in the graph of cycles and 0 if not
int is_edge(GRAPHE_CYCLE cy, int x, int y);


// Function to check if the given set of vertices
// in store array is a clique or not
int is_clique(GRAPHE_CYCLE cy, int b, int* store);

// Returns the address of the vertex with the given id, or NULL if does not exist
SOMMET* sommet_by_id(GRAPHE_CYCLE cy, int id);

// Returns 1 if the edges (s11, s12), (s21, s22) and (s31, s32) are identical
// (s11 and s12 are two consecutive vertices in one cycle, s21 and s22 two consecutive vertices in another cycle, and idem for s31 and s32) 
int meme_arete(int s11, int s12, int s21, int s22, int s31, int s32);

// Returns 1 if the edges (s11, s12) and  (s21, s22) are identical
// (s11 and s12 are two consecutive vertices in one cycle, s21 and s22 two consecutive vertices in another cycle) 
int meme_arete_2_cycles(int s11, int s12, int s21, int s22);

// Returns 1 if the vertex s is contained in the tabular tab_sommet of the size nb_sommet
// and 0 otherwise
int sommet_in_tab(int s, int* tab_sommet, int nb_sommets);

// Returns 1 if the edge (sommet1, sommet2) exists in the molecular graph m 
// and 0 otherwise
int arete_dans_molecule(struct molecule m, int sommet1, int sommet2);

// Returns the type (1: simple bond or 2: double bond) of the edge (sommet1, sommet2) if it exists in the molecular graph m
// and -1 otherwise
int type_arete_dans_molecule(struct molecule m, int sommet1, int sommet2);

// Returns the number of common vertices between two cliques in the graph of cycles (stored in the "ids" field of vertices in the graph of corners)
int nb_common_cycles(SOMMET_COIN sommet1, SOMMET_COIN sommet2, int taille_cliques);

// Returns 1 if the vertex with the id num_sommet is in a bouboule 
// and 0 otherwise 
int sommet_dans_bouboule(GRAPHE_CYCLE cy, int num_sommet, int* bouboule, int taille_bouboule);

