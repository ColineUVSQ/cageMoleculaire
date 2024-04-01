#ifndef NOM_D3
#define NOM_D3

#include "structure.h"

void init_atom_num ();
struct molecule lire_molecule_mol(FILE *F);
double chrono();
int lire_num_atome(FILE *F);
int valeur_char (FILE *F) ;
void ligne_suivante(FILE *F);
int lire_entier_3 (FILE * F);
struct liaison lire_liaison(FILE *F);
void lire_fin_molecule(FILE *F);
void trouver_la_fin_de_M(FILE *F);
void liberer_molecule(struct molecule m);
int atom_num (char *name);

#endif
