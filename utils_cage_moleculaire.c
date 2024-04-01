#include "utils_cage_moleculaire.h"

MOL_CARAC *classification = NULL;

// classification: name = nom de la molécule, classe = 1, 2, ou 3 pour cage, précage, non cage
void init_cage_non_cage()
{
	classification[0].name = "rugulin";
	classification[0].classe = 1;
	classification[1].name = "tsukushinamine_A";
	classification[1].classe = 2;
	classification[2].name = "artelein";
	classification[2].classe = 1;
	classification[3].name = "trichodimerol";
	classification[3].classe = 1;
	classification[4].name = "kopsinitarine_A";
	classification[4].classe = 1;
	classification[5].name = "sampsonione_A";
	classification[5].classe = 1;
	classification[6].name = "bielschowskysin";
	classification[6].classe = 2;
	classification[7].name = "scholarisine_A";
	classification[7].classe = 2;
	classification[8].name = "jiadifenolide";
	classification[8].classe = 2;
	classification[9].name = "platensimycin";
	classification[9].classe = 2;
	classification[10].name = "hippolachnin_A";
	classification[10].classe = 2;
	classification[11].name = "borneol";
	classification[11].classe = 3;
	classification[12].name = "atropurpuran";
	classification[12].classe = 3;
	classification[13].name = "jatamanin_H";
	classification[13].classe = 3;
	classification[14].name = "varioxepine_A";
	classification[14].classe = 3;
	classification[15].name = "illihenin_A";
	classification[15].classe = 3;
	classification[16].name = "arboduridine";
	classification[16].classe = 2;
	classification[17].name = "epicoccane_C";
	classification[17].classe = 2; // pre cage finalement
	classification[18].name = "kopsine";
	classification[18].classe = 3;
	classification[19].name = "harringtonolide";
	classification[19].classe = 2;
	classification[20].name = "daphniacetal_A";
	classification[20].classe = 2;
	classification[21].name = "alternarilactone_A";
	classification[21].classe = 2;
	classification[22].name = "nesteretal_A";
	classification[22].classe = 2;
	classification[23].name = "punctaporonin_U";
	classification[23].classe = 2;
	classification[24].name = "illisimonin_A";
	classification[24].classe = 2;
	classification[25].name = "astellolide_R";
	classification[25].classe = 3;
	classification[26].name = "arthpyrone_M";
	classification[26].classe = 3;
	classification[27].name = "aconicumine_A";
	classification[27].classe = 2;
	classification[28].name = "arsenicin_A";
	classification[28].classe = 1;
	classification[29].name = "euphylonoid_A";
	classification[29].classe = 2;
	classification[30].name = "bislangduoid_A";
	classification[30].classe = 2;
	classification[31].name = "alstoscholarisine_K";
	classification[31].classe = 1;
	classification[32].name = "granatripodin_B";
	classification[32].classe = 2;
	classification[33].name = "meloyinnanine_A";
	classification[33].classe = 2; // pre-cage finalement
	classification[34].name = "longifolactone_F";
	classification[34].classe = 2;
	classification[35].name = "aplydactone";
	classification[35].classe = 3;
	classification[36].name = "tabernabovine_B";
	classification[36].classe = 2;
	classification[37].name = "asperazine_D";
	classification[37].classe = 3;
	classification[38].name = "inaequalisine_A";
	classification[38].classe = 3;
	classification[39].name = "strepantibin_A";
	classification[39].classe = 3;
	classification[40].name = "KB343";
	classification[40].classe = 3;
	classification[41].name = "geopyxin_A";
	classification[41].classe = 3;
	classification[42].name = "lobarialide_A";
	classification[42].classe = 3;
	classification[43].name = "cycloaplysinopsin_A";
	classification[43].classe = 3;
	classification[44].name = "aleutianamine";
	classification[44].classe = 3;
	classification[45].name = "peshawaraquinone";
	classification[45].classe = 3;
	classification[46].name = "retigeranic_acid_A";
	classification[46].classe = 3;
	classification[47].name = "scabrolide_A";
	classification[47].classe = 2; // pre cage finalement
	classification[48].name = "greenwaylactam_A";
	classification[48].classe = 3;
	classification[49].name = "ineleganolide";
	classification[49].classe = 2;
	classification[50].name = "pedrolide";
	classification[50].classe = 2;
	classification[51].name = "anislactone_B";
	classification[51].classe = 3;
	classification[52].name = "ceforalide_F";
	classification[52].classe = 2;
	classification[53].name = "malibatol_A";
	classification[53].classe = 3;
	classification[54].name = "incargranine_A";
	classification[54].classe = 3;
	classification[55].name = "brevianamide_A";
	classification[55].classe = 3;
	classification[56].name = "ulodione_A";
	classification[56].classe = 3;
	classification[57].name = "stemoamide";
	classification[57].classe = 3;
	classification[58].name = "schweinfurthin_Q";
	classification[58].classe = 3;
	classification[59].name = "daphnodorin_A";
	classification[59].classe = 3;
	classification[60].name = "volvaltrate_A";
	classification[60].classe = 3;
	classification[61].name = "galanthamine";
	classification[61].classe = 3;
	classification[62].name = "havellockate";
	classification[62].classe = 2;
	classification[63].name = "myriberine";
	classification[63].classe = 3;
	classification[64].name = "myrioxazine_A";
	classification[64].classe = 3;
	classification[65].name = "epicoccin_G";
	classification[65].classe = 3;
	classification[66].name = "biyouyanagin_A";
	classification[66].classe = 3;
	classification[67].name = "aberrarone";
	classification[67].classe = 3;
	classification[68].name = "daphmanidin";
	classification[68].classe = 3;
	classification[69].name = "spirotryprostatin_B";
	classification[69].classe = 3;
	classification[70].name = "phorbol";
	classification[70].classe = 3;
	classification[71].name = "xiamycin_A";
	classification[71].classe = 3;
	classification[72].name = "okaramine_N";
	classification[72].classe = 3;
	classification[73].name = "chartelline_C";
	classification[73].classe = 3;
	classification[74].name = "vinigrol";
	classification[74].classe = 3;
	classification[75].name = "maoecrystal";
	classification[75].classe = 3;
	classification[76].name = "pallambin_C";
	classification[76].classe = 3;
	classification[77].name = "taxuyunnanine_D";
	classification[77].classe = 3;
	classification[78].name = "ingenol";
	classification[78].classe = 3;
	classification[79].name = "cortistatin_A";
	classification[79].classe = 3;
	classification[80].name = "santonin";
	classification[80].classe = 3;
}



// Returns a 2 dimensions tabular containing nb_lignes and nb_col with -1 everywhere
int** init_tab_deux_dimensions(int nb_lignes, int nb_col)
{
	int i,j;
	int** tab =  malloc(nb_lignes*sizeof(int*));
	if(tab == NULL)
	{
		printf("Allocation échouée de tab bouboule\n");
		exit(4);
	}

	for(i = 0; i < nb_lignes; ++i)
	{
		tab[i] = malloc(nb_col*sizeof(int));
		if(tab[i] == NULL)
		{
			printf("Allocation échouée de tab bouboule %d\n", i);
			exit(5);
		}
		for(j=0;j<nb_col;j++)
		{
			tab[i][j] = -1;
		}
	}
	return tab;
}


// Return the vertex a in which are copied the values of fields of the vertex b
SOMMET copier_sommet(SOMMET a, SOMMET b)
{
	a.type = b.type;
	a.poids = b.poids;
	a.id = b.id;
	a.poids_bouboule = b.poids_bouboule;
	return a;
}

// Returns a tab with all the degrees of the vertices in the graph of cycles cy
int *calcul_degre(GRAPHE_CYCLE cy)
{
	int i; 
	int num_max = 0;
	for(i = 0; i < cy.nb_sommets; i++)
	{
		if (cy.liste_sommets[i].id > num_max) num_max = cy.liste_sommets[i].id;
	}
	int *degre = malloc((num_max+1) * sizeof( int)); // le degré de tous les sommets du graphe (num max pour le cas où les id des sommets ne sont pas de 0 à n)
	
	for(i = 0; i < (num_max+1); i++)
		degre[i] = 0;
		
	//remplissage
	int pos1,pos2;
	for ( i = 0; i < cy.nb_aretes; i++)
	{
		pos1 = cy.liste_aretes[i].id1;
		pos2 = cy.liste_aretes[i].id2;
		degre[pos1]++;
		degre[pos2]++;
	}

	return degre;
}

// Function that returns 1 if (x,y) is an edge in the graph of cycles and 0 if not
int is_edge(GRAPHE_CYCLE cy, int x, int y)
{
	int i;
	for(i=0;i<cy.nb_aretes;i++)
	{
		if((cy.liste_aretes[i].id1 == x && cy.liste_aretes[i].id2 == y) || (cy.liste_aretes[i].id1 == y && cy.liste_aretes[i].id2 == x)){
			return 1;
		}
	}
	return 0; 
}

// Function to check if the given set of vertices
// in store array is a clique or not
int is_clique(GRAPHE_CYCLE cy, int b, int* store)
{
	int i,j;
    // Run a loop for all the set of edges
    // for the select vertex
    for (i = 0; i < b; i++) {
        for (j = i + 1; j < b; j++)
 
            // If any edge is missing
            if (is_edge(cy,store[i],store[j]) == 0)
                return 0;
    }
    return 1;
}

// Returns the address of the vertex with the given id, or NULL if does not exist
SOMMET* sommet_by_id(GRAPHE_CYCLE cy, int id)
{
	int i;
	for(i=0;i<cy.nb_sommets;i++)
	{
		if(cy.liste_sommets[i].id == id)
		{
			return &cy.liste_sommets[i];
		}
	}
	return NULL;
}

// Returns 1 if the edges (s11, s12), (s21, s22) and (s31, s32) are identical
// (s11 and s12 are two consecutive vertices in one cycle, s21 and s22 two consecutive vertices in another cycle, and idem for s31 and s32) 
int meme_arete(int s11, int s12, int s21, int s22, int s31, int s32)
{
	if (((s11 == s21 && s12 == s22) || (s11 == s22 && s12 == s21)) && ((s11 == s31 && s12 == s32) || (s11 == s32 && s12 == s31)))
	{
		return 1;
	} 
	return 0;
}  


// Returns 1 if the edges (s11, s12) and  (s21, s22) are identical
// (s11 and s12 are two consecutive vertices in one cycle, s21 and s22 two consecutive vertices in another cycle) 
int meme_arete_2_cycles(int s11, int s12, int s21, int s22)
{
	if ((s11 == s21 && s12 == s22) || (s11 == s22 && s12 == s21)) return 1;
	return 0;
}

// Returns 1 if the vertex s is contained in the tabular tab_sommet of the size nb_sommet
// and 0 otherwise
int sommet_in_tab(int s, int* tab_sommet, int nb_sommets)
{
	int i = 0;
	while(i<nb_sommets && tab_sommet[i] != s)
	{
		i++;
	}
	if(i == nb_sommets) return 0;
	return 1;
}

// Returns 1 if the edge (sommet1, sommet2) exists in the molecular graph m 
// and 0 otherwise
int arete_dans_molecule(struct molecule m, int sommet1, int sommet2)
{
	int i;
	for(i=0;i<m.nb_liaisons;i++)
	{
		if((m.liste_liaisons[i].A1 == sommet1 && m.liste_liaisons[i].A2 == sommet2) || (m.liste_liaisons[i].A2 == sommet1 && m.liste_liaisons[i].A1 == sommet2))
		{
			return 1;
		}
	}
	return 0;
}

// Returns the type (1: simple bond or 2: double bond) of the edge (sommet1, sommet2) if it exists in the molecular graph m
// and -1 otherwise
int type_arete_dans_molecule(struct molecule m, int sommet1, int sommet2)
{
	int i;
	for(i=0;i<m.nb_liaisons;i++)
	{
		if((m.liste_liaisons[i].A1 == sommet1 && m.liste_liaisons[i].A2 == sommet2) || (m.liste_liaisons[i].A2 == sommet1 && m.liste_liaisons[i].A1 == sommet2))
		{
			//printf("A1 %d A2 %d type %d\n", m.liste_liaisons[i].A1, m.liste_liaisons[i].A2, m.liste_liaisons[i].l_type);
			return m.liste_liaisons[i].l_type;
		}
	}
	return -1;
}

// Returns the number of common vertices between two cliques in the graph of cycles (stored in the "ids" field of vertices in the graph of corners)
int nb_common_cycles(SOMMET_COIN sommet1, SOMMET_COIN sommet2, int taille_cliques)
{
	int i,j;
	int nb = 0;
	for(i=0;i<taille_cliques;i++)
	{
		for(j=0;j<taille_cliques;j++)
		{
			if(sommet1.ids[i] == sommet2.ids[j]) nb ++;
		}
	}
	return nb;
}


// Returns 1 if the vertex with the id num_sommet is in a bouboule 
// and 0 otherwise 
int sommet_dans_bouboule(GRAPHE_CYCLE cy, int num_sommet, int* bouboule, int taille_bouboule)
{
	int i;
	for(i=0;i<taille_bouboule;i++)
	{
		if(num_sommet == bouboule[i])
		{
			return 1;
		}
	}
	return 0;
}
