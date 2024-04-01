#ifndef MON_H
#define MON_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define NB_ATOM_NAMES 119

char *atom_name[NB_ATOM_NAMES];


struct graphe{
	int nb_sommets;
	int nb_arete;
	int nb_connexe;
	struct lesommet *som;
	struct couplet *aretes;
	int **matrice_cycles_type;
	int **matrice_cycles_poids;
}graphe;


struct molecule{
	int nb_atomes;
	int nb_hydrogene;
	int nb_liaisons;
	int *liste_atomes; // liste des types d'atomes représentés par leur numéro dans le tableau des num d'atomes
	int **matrice_liaisons;
	struct liaison *liste_liaisons;
	struct graphe g;
	int g_def;
} molecule;

struct nom_at { char c1, c2; };
struct liaison { int A1, A2; int l_type; };

struct cycle{
	int poids;
	int *c;
}cycle;

struct isthme
{
	struct liaison l;
	int id_composant;
}isthme;
typedef struct isthme isthmes;
struct couple
{
	int a1;
	int a2;

}couple;
struct uncycle
{
	int nb_atomes;
	int sommet;
	struct liaison arete;
	int *chemin1;
	int *chemin2;
	int pere;
	int *sommets;
	int id_cycle;
}uncycle;

typedef struct uncycle cycles;

struct liste_voisins{
	int id_atome;
	int nb_voisins;
	int *id_voisins;
		
}liste_voisins;

struct couplet{
	int a1;
	int a2;
	int poids;
	int type; //1 s'il s'agit d'une liaison distance et 0 si les deux cycles partagent des liaison
}couplet;

struct lesommet{
	int id;
	int taille;
}lesommet;


struct type_arete{
	int type;
	int poids;
}type_arete;
struct arete_base
{
	int id1;
	int id2;
	int type;
	int poids;
}arete_base;
typedef struct arete_base ARETE;

struct unsommet
{
	int id;
	int poids;
	int type; // si contraction bouboule (0) ou non (1)
	char* poids_bouboule;

}unsommet;
typedef struct unsommet SOMMET;

struct graphe_cycle
{
	int nb_sommets;
	int nb_aretes;
	ARETE *liste_aretes;
	SOMMET *liste_sommets;
}graphe_cycle;


struct graphemoleculaire
{
	int nb_atomes;
	int nb_liaisons;
	int *liste_atomes;
	int *type_atomes;
	struct liaison *liste_liaisons;
	int **matrice_liaisons;
	int nb_connexe;
	int pere;//son numero de generation : 0 graphe de la molécule initial 1 : composantes connexe du grape initial 2: sousgraphes obentus en retirant les isthmes
}graphemoleculaire;

typedef struct graphemoleculaire graphemol;

typedef struct graphe_cycle GRAPHE_CYCLE;
int nb_arete_base;
int taille_base;
ARETE *base_aretes;
int *arete_cycle;
int **arete_liste;
cycles *labase;


// structure pour le graphe Ur
struct arc
{
	int id1;
	int id2;
};
typedef struct arc ARC;

struct sommet_vr
{
	int id;
};
typedef struct sommet_vr SOMMET_VR;

struct graphe_dr
{
	int nb_sommets;
	int nb_arcs;
	ARC *liste_arcs;
	SOMMET_VR *liste_sommets;
};


#endif
