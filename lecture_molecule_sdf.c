#include "lecture_molecule_sdf.h"

/**
 * Prétraitement des fichiers .mol ou .sdf
 * **/

void init_atom_num () {
	atom_name[  0] = "*";
	atom_name[  1] = "H";
	atom_name[  2] = "He";
	atom_name[  3] = "Li";
	atom_name[  4] = "Be";
	atom_name[  5] = "B";
	atom_name[  6] = "C";
	atom_name[  7] = "N";
	atom_name[  8] = "O";
	atom_name[  9] = "F";
	atom_name[ 10] = "Ne";
	atom_name[ 11] = "Na";
	atom_name[ 12] = "Mg";
	atom_name[ 13] = "Al";
	atom_name[ 14] = "Si";
	atom_name[ 15] = "P";
	atom_name[ 16] = "S";
	atom_name[ 17] = "Cl";
	atom_name[ 18] = "Ar";
	atom_name[ 19] = "K";
	atom_name[ 20] = "Ca";
	atom_name[ 21] = "Sc";
	atom_name[ 22] = "Ti";
	atom_name[ 23] = "V";
	atom_name[ 24] = "Cr";
	atom_name[ 25] = "Mn";
	atom_name[ 26] = "Fe";
	atom_name[ 27] = "Co";
	atom_name[ 28] = "Ni";
	atom_name[ 29] = "Cu";
	atom_name[ 30] = "Zn";
	atom_name[ 31] = "Ga";
	atom_name[ 32] = "Ge";
	atom_name[ 33] = "As";
	atom_name[ 34] = "Se";
	atom_name[ 35] = "Br";
	atom_name[ 36] = "Kr";
	atom_name[ 37] = "Rb";
	atom_name[ 38] = "Sr";
	atom_name[ 39] = "Y";
	atom_name[ 40] = "Zr";
	atom_name[ 41] = "Nb";
	atom_name[ 42] = "Mo";
	atom_name[ 43] = "Tc";
	atom_name[ 44] = "Ru";
	atom_name[ 45] = "Rh";
	atom_name[ 46] = "Pd";
	atom_name[ 47] = "Ag";
	atom_name[ 48] = "Cd";
	atom_name[ 49] = "In";
	atom_name[ 50] = "Sn";
	atom_name[ 51] = "Sb";
	atom_name[ 52] = "Te";
	atom_name[ 53] = "I";
	atom_name[ 54] = "Xe";
	atom_name[ 55] = "Cs";
	atom_name[ 56] = "Ba";
	atom_name[ 57] = "La";
	atom_name[ 58] = "Ce";
	atom_name[ 59] = "Pr";
	atom_name[ 60] = "Nd";
	atom_name[ 61] = "Pm";
	atom_name[ 62] = "Sm";
	atom_name[ 63] = "Eu";
	atom_name[ 64] = "Gd";
	atom_name[ 65] = "Tb";
	atom_name[ 66] = "Dy";
	atom_name[ 67] = "Ho";
	atom_name[ 68] = "Er";
	atom_name[ 69] = "Tm";
	atom_name[ 70] = "Yb";
	atom_name[ 71] = "Lu";
	atom_name[ 72] = "Hf";
	atom_name[ 73] = "Ta";
	atom_name[ 74] = "W";
	atom_name[ 75] = "Re";
	atom_name[ 76] = "Os";
	atom_name[ 77] = "Ir";
	atom_name[ 78] = "Pt";
	atom_name[ 79] = "Au";
	atom_name[ 80] = "Hg";
	atom_name[ 81] = "Tl";
	atom_name[ 82] = "Pb";
	atom_name[ 83] = "Bi";
	atom_name[ 84] = "Po";
	atom_name[ 85] = "At";
	atom_name[ 86] = "Rn";
	atom_name[ 87] = "Fr";   
	atom_name[ 88] = "Ra";
	atom_name[ 89] = "Ac";
	atom_name[ 90] = "Th";
	atom_name[ 91] = "Pa";
	atom_name[ 92] = "U";
	atom_name[ 93] = "Np";
	atom_name[ 94] = "Pu";
	atom_name[ 95] = "Am";
	atom_name[ 96] = "Cm";
	atom_name[ 97] = "Bk";
	atom_name[ 98] = "Cf";
	atom_name[ 99] = "Es";
	atom_name[100] = "Fm";
	atom_name[101] = "Md";
	atom_name[102] = "No";
	atom_name[103] = "Lr";
	atom_name[104] = "Rf";
	atom_name[105] = "Db";
	atom_name[106] = "Sg";
	atom_name[107] = "Bh";
	atom_name[108] = "Hs";
	atom_name[109] = "Mt";
	atom_name[110] = "Ds";
	atom_name[111] = "Rg";
	atom_name[112] = "Cn";
	atom_name[113] = "Uut";
	atom_name[114] = "Fl";
	atom_name[115] = "Uup";
	atom_name[116] = "Lv";
	atom_name[117] = "Uus";
	atom_name[118] = "Uuo";
}


int atom_num (char *name) {
	int i;
	for (i=0 ; i<NB_ATOM_NAMES ; i++) {
		if (strcmp(atom_name[i],name)==0){
			return i;}
		
	}
	return 0;
}

struct molecule lire_molecule_mol(FILE *F) 
{
	struct molecule M;
	// On saute l'entête
	ligne_suivante(F);
	ligne_suivante(F);
	ligne_suivante(F);
	if(fgetc(F) != '\n') // la 4e ligne n'est pas vide
	{
		fseek(F, -1, SEEK_CUR);
	}
	
	// Nombre d'atomes et nombre de liaisons
	M.nb_atomes = lire_entier_3(F);
	M.nb_liaisons = lire_entier_3(F);
	//printf("%d %d\n",M.nb_atomes,M.nb_liaisons);
	ligne_suivante(F);
	// Allocation mémoire
	M.liste_atomes   = malloc(M.nb_atomes   * sizeof(int));
	
	if(M.liste_atomes == NULL)
	{
		fprintf(stdout,"Cannot allocate memory  for M.liste_atomes\n");
		exit(55);
	}
	M.liste_liaisons = malloc(M.nb_liaisons * sizeof(struct liaison));
	
	if(M.liste_liaisons == NULL)
	{
		fprintf(stdout,"Cannot allocate memory  for M.liste_liaisons\n");
		exit(56);
	}
	int l,i,j; 
	M.matrice_liaisons = malloc(M.nb_atomes *sizeof(int *));
	for (l = 0; l < M.nb_atomes; l++)
		M.matrice_liaisons[l] = malloc( M.nb_atomes *sizeof(int));
	
	M.g.som = NULL;
	M.g.aretes = NULL;
	M.g.matrice_cycles_type = NULL;
	M.g.matrice_cycles_poids = NULL;
	M.nb_hydrogene=0;
	M.g_def = 0;
	M.g.nb_connexe = 0;
	// Lecture du nom des atomes
	int a; for (a=0 ; a<M.nb_atomes ; a++) 
	{
		
		M.liste_atomes[a] = lire_num_atome(F);
		//printf("%d\n", a);
		if(M.liste_atomes[a] == 1) M.nb_hydrogene++;
		ligne_suivante(F);
		
	}
	
	//printf("fin atomes\n");
	// Lecture des liaisons
	for (l=0 ; l<M.nb_liaisons ; l++)
	{
		M.liste_liaisons[l] = lire_liaison(F);
		if(M.liste_liaisons[l].A1 > M.nb_atomes || M.liste_liaisons[l].A2 > M.nb_atomes){ fprintf(stderr,"numero de molecule non valide \n"); exit(4); }
	}

	lire_fin_molecule(F);
	//remplissage de la matrice liaison
	
	for (i=0 ; i<M.nb_atomes ; i++)
	{
		for (j=0 ; j<M.nb_atomes; j++)
			M.matrice_liaisons[i][j] = 0;
	}
	for (l=0 ; l<M.nb_liaisons ; l++)
	{
		//printf("%d %d\n",M.liste_liaisons[l].A1,M.liste_liaisons[l].A2 );
		M.matrice_liaisons[M.liste_liaisons[l].A1 -1][M.liste_liaisons[l].A2 -1] = M.liste_liaisons[l].l_type;
		M.matrice_liaisons[M.liste_liaisons[l].A2 -1][M.liste_liaisons[l].A1 -1] = M.liste_liaisons[l].l_type;
	}

	return M;
}

double chrono() 
{
	struct timeval tv;
	static double date_deb = 0.0;
	double date_courante;
	gettimeofday(&tv,NULL);
	date_courante = tv.tv_sec + ((double)tv.tv_usec)*1e-6;
	if (date_deb == 0.0) date_deb = date_courante;
	return date_courante-date_deb;
}

int lire_num_atome(FILE *F) {
	// On saute les coordonnées
	float x,y,z;
	fscanf(F,"%f",&x);
	fscanf(F,"%f",&y);
	fscanf(F,"%f",&z);
	//printf("%f %f %f\n", x, y, z);
	// On saute un espace
	fgetc(F);
	// On lit le nom de l'atome
	char name[8];
	int c;
	int i = 0;
	c = fgetc(F);
	//printf("%c \n", c);
	while (c != ' ') {
		name[i] = c;
		i++;
		if (i>=8) {fprintf(stderr," pb dans lire_atome i>=8 : i vaut : %d\n",i); exit(4);}
		c= fgetc(F);
		//printf("%c\n", c);
	}
	
	// On termine la chaîne
	name[i] = '\0';
	//printf("%s\n", name);
	//printf("%d\n", atom_num(name));
	// On renvoie le numéro atomique
	return atom_num(name);
}


int valeur_char (FILE *F) {
	char c = fgetc(F);
	if (c==' ') return 0;
	if ((c<'0') || ('9'<c)) { fprintf(stderr,"pb caractère de nombre non valide : %d #%c#\n",c,c); exit(2); }
	return c-'0';
}

void ligne_suivante(FILE *F){
	static int nb_lignes = 0;
	while(fgetc(F) != '\n');
	nb_lignes++;
	
	//printf("nbl : %d\n",nb_lignes);
}

int lire_entier_3 (FILE * F) {
	char c1 = valeur_char(F);
	char c2 = valeur_char(F);
	char c3 = valeur_char(F);
	return 100*c1 + 10*c2 + c3;
}

struct liaison lire_liaison(FILE *F) {
	struct liaison l;
	l.A1 = lire_entier_3(F);
	l.A2 = lire_entier_3(F);
	l.l_type = lire_entier_3(F);
	//printf("\tliaison : %d %d %d\n",l.A1,l.A2,l.l_type);
	ligne_suivante(F);
	return l;
}


void lire_fin_molecule(FILE *F) {
	while (fgetc(F) != '$') ;
	ligne_suivante(F);
}

void trouver_la_fin_de_M(FILE *F) {
	while (fgetc(F)!='>')
		ligne_suivante(F);
}
void liberer_molecule(struct molecule m)
{

	if( m.liste_atomes != NULL)
		free(m.liste_atomes);
	if ( m.liste_liaisons != NULL)
		free(m.liste_liaisons);
	int i;
	if ( m.matrice_liaisons != NULL)
	{
		for ( i = 0; i < m.nb_atomes; i++)
			free(m.matrice_liaisons[i]);
		free(m.matrice_liaisons);
	}
	if(m.g.som != NULL)
		free(m.g.som);
	if(m.g.aretes != NULL)
		free(m.g.aretes);
	if(m.g.matrice_cycles_type != NULL)
	{
		for ( i = 0; i < m.g.nb_sommets; i++)
		{
			free(m.g.matrice_cycles_type[i]);
		}
		free(m.g.matrice_cycles_type);
	}

	if(m.g.matrice_cycles_poids != NULL)
	{
		for ( i = 0; i < m.g.nb_sommets; i++)
		{
			free(m.g.matrice_cycles_poids[i]);
		}
		free(m.g.matrice_cycles_poids);
	}

}


