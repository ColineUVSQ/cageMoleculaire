#include "analyse_cage.h"

/**
 * Calcul et étude des graphes de coins
 * 
 * */


// Returns the graph of cycles in which are removed the edges of type long and of weight 0
GRAPHE_CYCLE remove_edges(GRAPHE_CYCLE cy)
{
	int nb_garder = 0; // nb d'arêtes à garder (type 1 et poids > 0)
	int i;
	for ( i = 0; i < cy.nb_aretes; i++)
	{
		if(cy.liste_aretes[i].type == 1 && cy.liste_aretes[i].poids > 0)
		{
			nb_garder ++;
		}
	}
	//printf("%d\n", nb_garder);
	ARETE* new_liste_aretes = malloc(nb_garder * sizeof(ARETE));

	int j = 0;
	for ( i = 0; i < cy.nb_aretes; i++)
	{
		if(cy.liste_aretes[i].type == 1 && cy.liste_aretes[i].poids > 0)
		{
			//printf("i:%d\n", i);
			new_liste_aretes[j] = copier_arete(new_liste_aretes[j], cy.liste_aretes[i]);
			j++;
		}
	}
	
	free(cy.liste_aretes);
	cy.liste_aretes = new_liste_aretes;
	cy.nb_aretes = nb_garder;
	return cy;
}

// Returns 1 if the 3-sized clique stored in store (id des noeuds de la clique dans store[0], store[1], store[2]) is a "bouboule"
// and 0 otherwise
int is_bouboule(GRAPHE_CYCLE cy, int* store)
{
	int i;
	int poids_s1, poids_s2, poids_s3;
	for(i=0;i<cy.nb_sommets;i++)
	{
		if(cy.liste_sommets[i].id == store[0]) poids_s1 = cy.liste_sommets[i].poids;
		else if (cy.liste_sommets[i].id == store[1]) poids_s2 = cy.liste_sommets[i].poids;
		else if (cy.liste_sommets[i].id == store[2]) poids_s3 = cy.liste_sommets[i].poids;
	}
	
	int poids_a1, poids_a2, poids_a3;
	

	for(i=0;i<cy.nb_aretes;i++)
	{
		if((cy.liste_aretes[i].id1 == store[0] && cy.liste_aretes[i].id2 == store[1]) || (cy.liste_aretes[i].id1 == store[1] && cy.liste_aretes[i].id2 == store[0]))
		{
			poids_a1 = cy.liste_aretes[i].poids;
		}
		else if((cy.liste_aretes[i].id1 == store[0] && cy.liste_aretes[i].id2 == store[2]) || (cy.liste_aretes[i].id1 == store[2] && cy.liste_aretes[i].id2 == store[0]))
		{
			poids_a2 = cy.liste_aretes[i].poids;
		}
		else if((cy.liste_aretes[i].id1 == store[1] && cy.liste_aretes[i].id2 == store[2]) || (cy.liste_aretes[i].id1 == store[2] && cy.liste_aretes[i].id2 == store[1]))
		{
			poids_a3 = cy.liste_aretes[i].poids;
		}
	}
	
	// if the weight of each vertex is equal to the sum of the weight of its incident edges then the clique is a bouboule
	if(poids_s1 == poids_a1 + poids_a2 && poids_s2 == poids_a1 + poids_a3 && poids_s3 == poids_a2 + poids_a3)
	{
		return 1;
	}
	return 0;
}

// Dans le tableau sommets stockant les sommets appartenant à un cycle, les sommets peuvent ne pas être dans l'ordre de parcours du cycle
// Si ce n'est pas le cas, on les met dans cet ordre
void verif_suite_sommets(struct molecule m, int* sommets, int nb_sommets)
{
	int i,j,k,c,s,existe,refaire;
	refaire = 0;
	for(i=0;i<nb_sommets;i++)
	{	
		existe = 0;
		for(j=0;j<m.nb_liaisons;j++)
		{
			
			if((m.liste_liaisons[j].A1 == sommets[i] && m.liste_liaisons[j].A2 == sommets[(i+1)%nb_sommets]) || (m.liste_liaisons[j].A2 == sommets[i] && m.liste_liaisons[j].A1 == sommets[(i+1)%nb_sommets]))
			{
				existe = 1;
			}
		}
		if(!existe)
		{
			refaire = 1;
		}
	}
	if(refaire)
	{
		int* sommets_temp = malloc(nb_sommets*sizeof(int));
		for(i=0;i<nb_sommets;i++)
		{
			sommets_temp[i] = -1;
		}
		c = 0;
		i = 0;
		while(c != -1)
		{
			s = sommets[c]; // numéro du sommet stocké dans la case c
			sommets_temp[i] = s;
			c = -1;
			k = 0;
			while(k<nb_sommets && c == -1)
			{
				if (!sommet_in_tab(sommets[k], sommets_temp, nb_sommets)) // si on n'a pas mis le keme sommet dans sommets_temp
				{
					for(j=0;j<m.nb_liaisons;j++) // on cherche s'il y a une liaison entre le ceme et le keme sommet du cycle
					{
						if((m.liste_liaisons[j].A1 == s && m.liste_liaisons[j].A2 == sommets[k]) || (m.liste_liaisons[j].A2 == s && m.liste_liaisons[j].A1 == sommets[k]))
						{
							i++;
							c = k;
							break;
						}
					}
				}
				k++;
			}
		}
		// on remet la suite de sommets trouvés dans la variable sommets 
		for(i=0;i<nb_sommets;i++)
		{
			sommets[i] = sommets_temp[i];
		}

		free(sommets_temp);
	}	
}



// Function to print the cliques in csv files
// in the file f_out : the id and weight of the vertices of the clique
// in the file f_out_type : the type of the clique
void print(int n, int *store, GRAPHE_CYCLE cy, int type, FILE* f_out, FILE* f_out_type)
{
	int i,j,k;
	
	fprintf(f_out_type, "%d ",type); 
	
    for (i = 0; i < n; i++)
	{	
		j = 0;
		while(j<cy.nb_sommets && cy.liste_sommets[j].id != store[i])
		{
			j++;
		}
		if(j < cy.nb_sommets)
		{
			fprintf(f_out, "%d  ", cy.liste_sommets[j].poids);
			//printf("%d ",  cy.liste_sommets[j].id);
		}
    }
        
    for (i = 0; i < n; i++)
    {
		for (j = i+1; j < n; j++)
		{
			for (k = 0; k < cy.nb_aretes; k++)
			{
				if(cy.liste_aretes[k].id1 == store[i] && cy.liste_aretes[k].id2 == store[j]) 
					fprintf(f_out, "(%d:%d): %d ", store[i], store[j], cy.liste_aretes[k].poids);
				if(cy.liste_aretes[k].id1 == store[j] && cy.liste_aretes[k].id2 == store[i]) 
					fprintf(f_out, "(%d:%d): %d ", store[j], store[i], cy.liste_aretes[k].poids);
			}
			
		}
	}
	
    fprintf(f_out, ", ");
    fprintf(f_out_type, ", ");
    //printf("\n");
    
}



// Returns 1 if the 3-sized clique stored in store (id des noeuds de la clique dans store[0], store[1], store[2]) is a "coin bouboule" (presque artificiel)
// which means that the clique contains two cycles of a bouboule and one other cycle that is not in a bouboule
// and 0 otherwise
int is_coin_bouboule(int* store, int** tab_bouboule, int nb_bouboules)
{
	int i,j;
	int dans_bouboule;
	for(j=0;j<nb_bouboules;j++)
	{
		dans_bouboule = 0;
		
		for(i=0;i<3;i++)
		{
			if(store[i] == tab_bouboule[j][0] || store[i] == tab_bouboule[j][1]|| store[i] == tab_bouboule[j][2]) dans_bouboule++;
		}
		if (dans_bouboule == 2)
		{
			return 1;
		}
	}
	return 0;
}


// Returns 1 if the 3-sized clique stored in store (id des noeuds de la clique dans store[0], store[1], store[2]) is a "coin pas ouvert"
// which means that the 3 cycles have at least one common edge.
// Returns 0 otherwise
int coin_pas_ouvert(GRAPHE_CYCLE cy, int* store, struct molecule m) 
{
	int i,j,k;
	// on cherche les cycles dans la base de cycles
	//printf("taille base %d\n", taille_base);
	//printf("la base : %p\n", labase[0]);

	cycles* cycle1 = malloc(sizeof(cycles));
	cycles* cycle2 = malloc(sizeof(cycles));
	cycles* cycle3 = malloc(sizeof(cycles));
	
	for(i=0;i<taille_base;i++)
	{
		if(labase[i].id_cycle == store[0]) (*cycle1) = labase[i];
		else if(labase[i].id_cycle == store[1]) (*cycle2) = labase[i];
		else if(labase[i].id_cycle == store[2]) (*cycle3) = labase[i];
	}
	
	// met les sommets de cycle1, cycle2, cycle3 dans un ordre de parcours des cycles si pas déjà le cas
	verif_suite_sommets(m, cycle1->sommets, cycle1->nb_atomes);
	verif_suite_sommets(m, cycle2->sommets, cycle2->nb_atomes);
	verif_suite_sommets(m, cycle3->sommets, cycle3->nb_atomes);
	
	// on cherche les aretes en commun aux trois cycles
	int commun;
	commun = 0;
	for(i=0;i<cycle1->nb_atomes;i++)
	{
		if (!arete_dans_molecule(m, cycle1->sommets[i], cycle1->sommets[(i+1)%cycle1->nb_atomes])) // verifie que l'arete (i, i+1) existe bien (normalement oui avec la verif d'avant)
		{
			printf("problème: pas d'arête %d,%d", cycle1->sommets[i], cycle1->sommets[(i+1)%cycle1->nb_atomes]);
			exit(1);
		}
		for(j=0;j<cycle2->nb_atomes;j++)
		{
			if (!arete_dans_molecule(m, cycle2->sommets[j], cycle2->sommets[(j+1)%cycle2->nb_atomes])) 
			{
				printf("problème: pas d'arête %d,%d", cycle2->sommets[j], cycle2->sommets[(j+1)%cycle2->nb_atomes]);
				exit(1);
			}
			for(k=0;k<cycle3->nb_atomes;k++)
			{
				if (!arete_dans_molecule(m, cycle3->sommets[k], cycle3->sommets[(k+1)%cycle3->nb_atomes])) 
				{
					printf("problème: pas d'arête %d,%d", cycle3->sommets[k], cycle3->sommets[(k+1)%cycle3->nb_atomes]);
					exit(1);
				}
				if(meme_arete(cycle1->sommets[i], cycle1->sommets[(i+1)%cycle1->nb_atomes], cycle2->sommets[j], cycle2->sommets[(j+1)%cycle2->nb_atomes], cycle3->sommets[k], cycle3->sommets[(k+1)%cycle3->nb_atomes]))
				{
					commun ++;
				}
				 
			}
		}
	}
	
	free(cycle1);
	free(cycle2);
	free(cycle3);
	
	if(commun > 0)
	{
		return 1;
	}
	return 0;
}


// Returns the number of double bonds located at the junction of 2 cycles in the 3-sized clique stored in store (id des noeuds de la clique dans store[0], store[1], store[2])
int* nb_double_liaisons(GRAPHE_CYCLE cy, int* store, struct molecule m) 
{
	int i,j;
	int nb_dl_seul, nb_dl_commun;
	nb_dl_seul = 0;
	nb_dl_commun = 0;
	// on cherche les cycles dans la base de cycles
	//printf("taille base %d\n", taille_base);
	//printf("la base : %p\n", labase[0]);

	cycles* cycle1 = malloc(sizeof(cycles));
	cycles* cycle2 = malloc(sizeof(cycles));
	cycles* cycle3 = malloc(sizeof(cycles));
	
	for(i=0;i<taille_base;i++)
	{
		if(labase[i].id_cycle == store[0]) (*cycle1) = labase[i];
		else if(labase[i].id_cycle == store[1]) (*cycle2) = labase[i];
		else if(labase[i].id_cycle == store[2]) (*cycle3) = labase[i];
	}
	
	// met les sommets de cycle1, cycle2, cycle3 dans un ordre de parcours des cycles si pas déjà le cas
	verif_suite_sommets(m, cycle1->sommets, cycle1->nb_atomes);
	verif_suite_sommets(m, cycle2->sommets, cycle2->nb_atomes);
	verif_suite_sommets(m, cycle3->sommets, cycle3->nb_atomes);
	
	int trouve = 0;
	// on cherche les double liaisons (un peu long, à arranger)
	for(i=0;i<cycle1->nb_atomes;i++)
	{
		if(type_arete_dans_molecule(m, cycle1->sommets[i], cycle1->sommets[(i+1)%cycle1->nb_atomes]) == 2) // double liaison
		{
			trouve = 0;
			j = 0;
			while(j<cycle2->nb_atomes && !trouve)
			{
				if(meme_arete_2_cycles(cycle1->sommets[i], cycle1->sommets[(i+1)%cycle1->nb_atomes], cycle2->sommets[j], cycle2->sommets[(j+1)%cycle2->nb_atomes]))
				{
					nb_dl_commun ++;
					trouve = 1;
				}
				j++;
			}
			j = 0;
			while(j<cycle3->nb_atomes && !trouve)
			{
				if(meme_arete_2_cycles(cycle1->sommets[i], cycle1->sommets[(i+1)%cycle1->nb_atomes], cycle3->sommets[j], cycle3->sommets[(j+1)%cycle3->nb_atomes]))
				{
					nb_dl_commun ++;
					trouve = 1;
				}
				j++;
			}
			if (!trouve) nb_dl_seul ++;
		}
	}
	
	for(i=0;i<cycle2->nb_atomes;i++)
	{
		if(type_arete_dans_molecule(m, cycle2->sommets[i], cycle2->sommets[(i+1)%cycle2->nb_atomes]) == 2) // double liaison
		{
			trouve = 0;
			j = 0;
			while(j<cycle1->nb_atomes && !trouve)
			{
				if(meme_arete_2_cycles(cycle2->sommets[i], cycle2->sommets[(i+1)%cycle2->nb_atomes], cycle1->sommets[j], cycle1->sommets[(j+1)%cycle1->nb_atomes]))
				{
					trouve = 1;
				}
				j++;
			}
			j = 0;
			while(j<cycle3->nb_atomes && !trouve)
			{
				if(meme_arete_2_cycles(cycle2->sommets[i], cycle2->sommets[(i+1)%cycle2->nb_atomes], cycle3->sommets[j], cycle3->sommets[(j+1)%cycle3->nb_atomes]))
				{
					nb_dl_commun ++;
					trouve = 1;
				}
				j++;
			}
			if (!trouve) nb_dl_seul ++;
		}
	}
	
	for(i=0;i<cycle3->nb_atomes;i++)
	{
		if(type_arete_dans_molecule(m, cycle3->sommets[i], cycle3->sommets[(i+1)%cycle3->nb_atomes]) == 2) // double liaison
		{
			trouve = 0;
			j = 0;
			while(j<cycle1->nb_atomes && !trouve)
			{
				if(meme_arete_2_cycles(cycle3->sommets[i], cycle3->sommets[(i+1)%cycle3->nb_atomes], cycle1->sommets[j], cycle1->sommets[(j+1)%cycle1->nb_atomes]))
				{
					trouve = 1;
				}
				j++;
			}
			j = 0;
			while(j<cycle2->nb_atomes && !trouve)
			{
				if(meme_arete_2_cycles(cycle3->sommets[i], cycle3->sommets[(i+1)%cycle3->nb_atomes], cycle2->sommets[j], cycle2->sommets[(j+1)%cycle2->nb_atomes]))
				{
					trouve = 1;
				}
				j++;
			}
			if (!trouve) nb_dl_seul ++;
		}
	}
	
	
	free(cycle1);
	free(cycle2);
	free(cycle3);
	
	int* nb_dl = malloc(2*sizeof(int));
	nb_dl[0] = nb_dl_commun; // commun à 2 cycles
	nb_dl[1] = nb_dl_seul; // appartient à un seul cycle
	return nb_dl;
	
}


// Find all the cliques of size s in the graph of cycles cy, and find bouboules among them and store them in tab_bouboule
/* i : starting node
 * l : length of the present set of nodes
 * s : size of clique
 * degre : all the degrees of the vertices of the graph
 * store : contains the current set of nodes that could form a clique
 * tab_bouboule : tableau à deux dimensions pour stocker les sommets des bouboules
 * b : nb of bouboules found so far
	Returns the number of bouboules found
 */
int findCliquesToDetectBouboules(GRAPHE_CYCLE cy, int i, int l, int s, int *degre, int *store,  int** tab_bouboule, int b)
{
	int n = cy.nb_sommets;
	int j;
    for (j = i; j < n - (s - l); j++)
	{ 
        if (degre[cy.liste_sommets[j].id] >= s - 1) {
 
            // Add the vertex to store
            store[l-1] = cy.liste_sommets[j].id;

            // check if the present set of nodes with node j can be a clique
            if (is_clique(cy, l, store)){
 
                // If the length of the clique is
                // still less than the desired size
                if (l < s){
 
                    // Recursion to add vertices
                    b = findCliquesToDetectBouboules(cy, j+1, l + 1, s, degre, store, tab_bouboule, b);
                }
                // Size is met
                else{
					if (is_bouboule(cy, store)){
						tab_bouboule[b][0] = store[0];
						tab_bouboule[b][1] = store[1];
						tab_bouboule[b][2] = store[2];
						b ++;
					}   
				}
            }
        }
     }
     return b;
}


// Find all the cliques of size s in the graph of cycles cy, and classify them in types and write the types in csv files
/* i : starting node
 * l : length of the present set of nodes
 * s : size of clique
 * degre : all the degrees of vertices of cy
 * store : contains the current set of nodes that could form a clique
 * f_out, f_out_type, f_out_dl : csv files (f_out to store informations about each clique (weight of vertices and edges for example), f_out_type to store the types of every clique, f_out_dl to store the number of double bonds at the junctions of cycles or not at the junction in each clique)
 * tab_bouboule : 2 dimensions tabular in which are stored the bouboules
 * b : number of bouboules found
 * c : current number of cliques of type 3 (not bouboule, not coin bouboule, not coin pas ouvert, no double bonds at the junction of cycles), 
 * m : the molecular graph
 * Returns the number of cliques we will put in the graph of corners (c)
 */
int findCliquesToCount(GRAPHE_CYCLE cy, int i, int l, int s, int *degre, int *store, FILE* f_out, FILE* f_out_type, FILE* f_out_dl, int** tab_bouboule, int b, int c, struct molecule m)
{
	int n = cy.nb_sommets;
	int j;
	int type;
	int* nb_dl;
    for (j = i; j < n - (s - l); j++)
	{ 
        if (degre[cy.liste_sommets[j].id] >= s - 1) {
 
            // Add the vertex to store
            store[l-1] = cy.liste_sommets[j].id;

            // check if the present set of nodes with node j can be a clique
            if (is_clique(cy, l, store)){
 
                // If the length of the clique is
                // still less than the desired size
                if (l < s){
 
                    // Recursion to add vertices
                    c = findCliquesToCount(cy, j+1, l + 1, s, degre, store, f_out, f_out_type, f_out_dl, tab_bouboule, b, c, m);
                }
                // Size is met
                else{
					
						if (is_bouboule(cy, store))
						{
							type = 0;
							//printf("type 0 : %d %d %d\n", store[0], store[1], store[2]);
						}
						else if (is_coin_bouboule(store, tab_bouboule, b))
						{
							type = 1;
							//printf("type 1 : %d %d %d\n", store[0], store[1], store[2]);
						}
						else if(coin_pas_ouvert(cy, store, m))
						{
							type = 2;
							//printf("type 2 : %d %d %d\n", store[0], store[1], store[2]);
						}
						else
						{
							nb_dl = nb_double_liaisons(cy, store, m);
							fprintf(f_out_dl, "(%d: %d),", nb_dl[0], nb_dl[1]);
							if(nb_dl[0] > 0)
							{
								type = 4; // des double liaisons aux jonctions
								//printf("type 4 : %d %d %d\n", store[0], store[1], store[2]);
							}
							else 
							{
								type = 3;
								//printf("type 3 : %d %d %d\n", store[0], store[1], store[2]);
								c++;
							}
							free(nb_dl);
						}
						print(l, store, cy, type, f_out, f_out_type);                    
				}
            }
        }
     }
     return c;
}


// Find all the cliques of size s and construct vertices of graph corner
/* i : starting node
 * l : length of the present set of nodes
 * s : size of clique
 * degre : all the degrees of vertices of cy
 * store : contains the current set of nodes that could form a clique
 * tab_bouboule : 2 dimensions tabular in which are stored the bouboules
 * b : number of bouboules found
 * c : current number of cliques of type 3 (not bouboule, not coin bouboule, not coin pas ouvert, no double bonds at the junction of cycles) = number of vertices in our graph of corners
 * liste_sommets : the list of vertices in progress that will contain the graph of corners
 * mol : the molecular graph
 * Returns the number of cliques we will put in the graph of corners (c)
 */
int findCliquesToConstructCoinGraph(GRAPHE_CYCLE cy, int i, int l, int s, int *degre, int *store, int** tab_bouboule, int b, int c, SOMMET_COIN* liste_sommets, struct molecule mol)
{
	int n = cy.nb_sommets;
	int j,k,m,indice;
	int* nb_dl;
    // Check if any vertices from i can be inserted
    for (j = i; j < n - (s - l); j++)
	{ 
        // If the degree of the node j is sufficient 
        //printf("j: %d\n", cy.liste_sommets[j].id);
        if (degre[cy.liste_sommets[j].id] >= s - 1) {
 
            // Add the vertex to store
            store[l-1] = cy.liste_sommets[j].id;

            // check if the present set of nodes with node j can be a clique
            if (is_clique(cy, l, store)){
 
                // If the length of the clique is
                // still less than the desired size
                if (l < s){
 
                    // Recursion to add vertices
                    c = findCliquesToConstructCoinGraph(cy, j+1, l + 1, s, degre, store, tab_bouboule, b, c, liste_sommets, mol);
                }
                // Size is met
                else{
						
						if (!is_bouboule(cy,store) && !is_coin_bouboule(store, tab_bouboule, b) && !coin_pas_ouvert(cy, store, mol)) 
						{
							nb_dl = nb_double_liaisons(cy, store, mol);
							if(nb_dl[0] == 0)
							{
								// c'est un potentiel coin de cage (type 3)
								liste_sommets[c].id = c;
								liste_sommets[c].ids = malloc(s*sizeof(int));
								liste_sommets[c].poids = malloc(s*sizeof(int));
								liste_sommets[c].liaisons_communs = malloc(s*sizeof(int));
								for(k=0;k<s;k++)
								{
									liste_sommets[c].ids[k] = store[k];
									indice = -1;
									for(m=0;m<cy.nb_sommets;m++)
									{
										if(store[k] == cy.liste_sommets[m].id)
										{
											indice = m;
										}
									}
									liste_sommets[c].poids[k] = cy.liste_sommets[indice].poids;
									
									// cherche le nb de liaisons entre les cycles de la clique
									for(m=0;m < cy.nb_aretes;m++)
									{
										if((cy.liste_aretes[m].id1 == store[(k+1)%s] && cy.liste_aretes[m].id2 == store[(k+2)%s]) || (cy.liste_aretes[m].id2 == store[(k+1)%s] && cy.liste_aretes[m].id1 == store[(k+2)%s]) ) 
										{
											liste_sommets[c].liaisons_communs[k] = cy.liste_aretes[m].poids;
										}
									}
								}
								
								c++;
							}
							free(nb_dl);
						}
						    
				}
            }
        }
     }
     return c;
}



// Returns the graph of corners obtained from the graph of cycles cy
/* nb_sommets : number of vertices in the graph of corners
 * degre : all the degrees of vertices in the graph of cycles
 * store : in the search for cliques, to store sets of nodes that could form a clique
 * tab_bouboule : 2 dimensions tabular in which are stored the bouboules
 * b : number of bouboules found
 * m : the molecular graph
 * */
GRAPHE_COIN construct_graph_coin(GRAPHE_CYCLE cy, int nb_sommets, int* degre, int* store, int** tab_bouboule, int b, struct molecule m)
{
	GRAPHE_COIN gc;
	int i,j,k;
	// creation des sommets
	gc.nb_sommets = nb_sommets;
	gc.liste_sommets = malloc(nb_sommets * sizeof(SOMMET_COIN));
	// pour remplir gc.liste_sommets : un sommet = un coin de type 3 (coin ouvert sans double liaisons aux jonctions des cycles)
	findCliquesToConstructCoinGraph(cy, 0, 1, 3, degre, store, tab_bouboule, b, 0, gc.liste_sommets, m);
	
	// pour remplir gc.liste_aretes : une arête entre deux coins s'ils ont un cycle en commun (pondérées par le nombre de cycles en commun)
	int nb_aretes = 0;
	for(i=0;i<gc.nb_sommets;i++)
	{
		for(j=i+1;j<gc.nb_sommets;j++)
		{
			if (nb_common_cycles(gc.liste_sommets[i], gc.liste_sommets[j], 3) > 0) nb_aretes ++;
		}
	}

	gc.liste_aretes = malloc(nb_aretes*(sizeof(ARETE_COIN)));
	gc.nb_aretes = nb_aretes;
	int poids_aretes;
	k = 0;
	for(i=0;i<gc.nb_sommets;i++)
	{
		for(j=i+1;j<gc.nb_sommets;j++)
		{
			poids_aretes = nb_common_cycles(gc.liste_sommets[i], gc.liste_sommets[j], 3);
			if(poids_aretes > 0)
			{
				gc.liste_aretes[k].id1 = i;
				gc.liste_aretes[k].id2 = j;
				gc.liste_aretes[k].poids = poids_aretes;
				k++;
			}
		}
	}
	return gc;
}



// Print the informations of a molecular graph
void afficheInfosGrapheMol(struct molecule m)
{
	int j;
	printf("nb atomes : %d, nb liaisons: %d\n", m.nb_atomes, m.nb_liaisons);
	for(j=0;j<m.nb_liaisons;j++)
	{
		printf("arete %d: id1 %d id2 %d type %d \n", j, m.liste_liaisons[j].A1, m.liste_liaisons[j].A2,  m.liste_liaisons[j].l_type); 
	}
	for(j=0;j<m.nb_atomes;j++)
	{
		printf("sommet %d: id %d \n", j, m.liste_atomes[j]); 
	}
}


// Print the informations of a graph of cycles
void afficheInfosGrapheCycle(GRAPHE_CYCLE cy)
{
	int j;
	printf("nb sommets : %d, nb aretes: %d\n", cy.nb_sommets, cy.nb_aretes);
	for(j=0;j<cy.nb_aretes;j++)
	{
		printf("arete %d: id1 %d id2 %d type %d poids %d\n", j, cy.liste_aretes[j].id1, cy.liste_aretes[j].id2,  cy.liste_aretes[j].type,  cy.liste_aretes[j].poids); 
	}
	for(j=0;j<cy.nb_sommets;j++)
	{
		printf("sommet %d: id %d type %d poids %d poids_bouboule %s\n", j, cy.liste_sommets[j].id, cy.liste_sommets[j].type,  cy.liste_sommets[j].poids,  cy.liste_sommets[j].poids_bouboule); 
	}
}

// Print the informations of a graph of corners
void afficheInfosGrapheCoin(GRAPHE_COIN c)
{
	int j,k;
	printf("Graphe de coin: \n nb sommets : %d, nb aretes: %d\n", c.nb_sommets, c.nb_aretes);
	for(j=0;j<c.nb_aretes;j++)
	{
		printf("arete %d: id1 %d id2 %d poids %d \n", j, c.liste_aretes[j].id1, c.liste_aretes[j].id2, c.liste_aretes[j].poids); 
	}
	for(j=0;j<c.nb_sommets;j++)
	{
		printf("sommet %d: id %d", j, c.liste_sommets[j].id);
		printf(" poids : ");
		for(k=0;k<3;k++)
		{
			printf("%d ",  c.liste_sommets[j].poids[k]); 
		}
		printf(" ids : ");
		for(k=0;k<3;k++)
		{
			printf("%d ",  c.liste_sommets[j].ids[k]); 
		}
		printf(" liaisons communes : ");
		for(k=0;k<3;k++)
		{
			printf("%d ",  c.liste_sommets[j].liaisons_communs[k]); 
		}
		printf("\n");
	}
}

void genere_dot_file_cycles(GRAPHE_CYCLE cy, char* name, int taille)
{
	int i;
	FILE* F1, *F2;
	int taille_file_dot = strlen(FILES_DOT_GCYCLES);
	char* name_file_dot = malloc((taille_file_dot+taille+1) * sizeof(char));
	strcat(strcpy(name_file_dot, FILES_DOT_GCYCLES), name);
	strcat(name_file_dot, ".dot");
	
	F1 = fopen(name_file_dot, "w");
	// Exporter pour Graphviz
	fprintf(F1, "graph G {\n");
	for(i=0;i<cy.nb_sommets;i++)
	{
		fprintf(F1, "\t%d [label=\"%d\"];\n", cy.liste_sommets[i].id, cy.liste_sommets[i].poids);
	} 
	for(i=0;i<cy.nb_aretes;i++)
	{	
		if (cy.liste_aretes[i].type == 1)
			fprintf(F1, "\t%d -- %d [label=\"%d\", color=blue];\n", cy.liste_aretes[i].id1, cy.liste_aretes[i].id2, cy.liste_aretes[i].poids);
		else 
			fprintf(F1, "\t%d -- %d [label=\"%d\", color=green];\n", cy.liste_aretes[i].id1, cy.liste_aretes[i].id2, cy.liste_aretes[i].poids);
	}
	fprintf(F1, "}");
	fclose(F1);
	
	// Exporter en png
	F2 = fopen("scripts/genere_images_g_cycles.sh","a");
	char* name_file_png = malloc((taille_file_dot+taille+1)*sizeof(char));
	strcat(strcpy(name_file_png, FILES_PNG_GCYCLES), name);
	strcat(name_file_png, ".png");
	
	fprintf(F2, "dot -Tpng %s -o %s\n", name_file_dot, name_file_png);
	fclose(F2);
	free(name_file_dot);
	free(name_file_png);
}

void genere_dot_file_coins(GRAPHE_COIN gc, char* name, int taille)
{
	FILE* F1, *F2;
	// Exporter pour Graphviz
	int taille_file_dot = strlen(FILES_DOT_GCOINS);
	char* name_file_dot = malloc((taille_file_dot+taille+1) * sizeof(char));
	strcat(strcpy(name_file_dot, FILES_DOT_GCOINS), name);
	strcat(name_file_dot, ".dot");
	
	int i,k;
	F1 = fopen(name_file_dot, "w");
	fprintf(F1, "graph G {\n");
	for(i=0;i<gc.nb_sommets;i++)
	{
		fprintf(F1, "\t%d [label=<", gc.liste_sommets[i].id);
		for(k=0;k<3;k++)
		{
			fprintf(F1, "%d ", gc.liste_sommets[i].poids[k]);
		}
		fprintf(F1, "<BR/> <FONT COLOR=\"BLUE\">");
		for(k=0;k<3;k++)
		{
			fprintf(F1, "%d ", gc.liste_sommets[i].liaisons_communs[k]);
		}
		fprintf(F1, "</FONT>>];\n");
		
	} 
	for(i=0;i<gc.nb_aretes;i++)
	{	
		fprintf(F1, "\t%d -- %d [label=\"%d\"];\n", gc.liste_aretes[i].id1, gc.liste_aretes[i].id2, gc.liste_aretes[i].poids);
	}
	fprintf(F1, "}");
	fclose(F1);
	F2 = fopen("scripts/genere_images_g_coins.sh","a");
	char* name_file_png = malloc((taille_file_dot+taille+1)*sizeof(char));
	strcat(strcpy(name_file_png, FILES_PNG_GCOINS), name);
	strcat(name_file_png, ".png");
	
	fprintf(F2, "dot -Tpng %s -o %s\n", name_file_dot, name_file_png);
	fclose(F2);
	free(name_file_png);
	free(name_file_dot);
}

// Free memory of a graph of corners
void liberer_graphe_coin(GRAPHE_COIN c)
{
	int i;
	if(c.liste_sommets != NULL){
		for(i=0;i<c.nb_sommets;i++)
		{
			free(c.liste_sommets[i].ids);
			free(c.liste_sommets[i].poids);
			free(c.liste_sommets[i].liaisons_communs);
		}
		free(c.liste_sommets);
	}
	if(c.liste_aretes != NULL)
		free(c.liste_aretes);
	
}

int main(int argc, char *argv[])
{
	struct dirent *lecture;
    DIR *rep;
    int i;
    
    char *name = NULL;
    char *name_file_mol = NULL;
	int taille;
    
    int taille_mol = strlen(FILES_SMI);
    
    /*
     * Fichiers csv pour les résultats
     * */
    FILE *F, *F_out, *F_out_type, *f_out_dl, *f_liste;
    F_out = fopen(RESULTS_ENS_CLIQUES, "w"); //liste des cliques par molécule avec les infos (poids des sommets et des arêtes)
    if(F_out == NULL)
    {
		printf("Impossible d'ouvrir le fichier %s\n", RESULTS_ENS_CLIQUES);
		exit(3);
	}
	F_out_type = fopen(RESULTS_TYPE_CLIQUES, "w"); // liste des cliques de chaque type par molécule
	if(F_out_type == NULL)
	{
		printf("Impossible d'ouvrir le fichier %s.csv\n", RESULTS_TYPE_CLIQUES);
		exit(4);
	}
	f_out_dl = fopen(RESULTS_DL_CLIQUES, "w"); //liste du nombre de double liaisons dans les cliques (2 valeurs: double liaisons communes à au moins 2 cycles, double liaisons dans un seul cycle)
	if(f_out_dl == NULL)
	{
		printf("Impossible d'ouvrir le fichier %s.csv\n", RESULTS_DL_CLIQUES);
		exit(5);
	}
	
	f_liste = fopen(RESULTS_LISTE_COINS, "w"); // liste des coins par molécule (poids des sommets et poids des arêtes)
	if(f_liste == NULL)
	{
		printf("Impossible d'ouvrir le fichier %s\n", RESULTS_LISTE_COINS);
		exit(6);
	}
	
	
    rep = opendir("data/smi_files_reduit"); // les .smi (SMILES notations) et les .mol (fichiers 3D pour input des constructions de graphes moléculaires) sont stockés dans le dossier smi_files_reduit
    
    // initialise la classification des chimistes en cage, precage, non cage
    classification = malloc(NB_MOL * sizeof(MOL_CARAC));
    init_cage_non_cage(); 
    
    while ((lecture = readdir(rep))) {
        printf("%s\n", lecture->d_name);
        // problème avec l'arsenichin A, SMILES non géré par OpenBabel
        if (strstr(lecture->d_name, ".mol") && strcmp(lecture->d_name, "arsenicin_A.mol")){ // !strcmp(lecture->d_name, "tabernabovine_B.mol")){
		
			taille = strlen(lecture->d_name);
			//printf("%d\n", taille);
			// allocation de memoire
			name = malloc((taille-3) * sizeof(char));
			name_file_mol = malloc((taille_mol+taille+1) * sizeof(char));
			
			// affectation des noms des fichiers
			strncpy(name, lecture->d_name, taille-4);
			name[taille-4] = '\0'; // il faut ajouter le caractère de fin de chaines

			strcat(strcpy(name_file_mol, FILES_SMI), lecture->d_name);
			
			// remplissage fichiers avec nom molécule et cage/précage/non cage
			F = fopen(name_file_mol,"r");
			if(F == NULL){
				printf("Impossible d'ouvrir le fichier %s", name_file_mol);
				exit(2);
			}
			
			fprintf(F_out, "%s,", name);
			fprintf(F_out_type, "%s,", name);
			fprintf(f_out_dl, "%s,", name);
			fprintf(f_liste, "%s,",name);
			
			for(i=0;i<NB_MOL;i++)
			{
				if(!strcmp(name, classification[i].name))
				{
					if(classification[i].classe == 1)
					{
						fprintf(F_out, "%s,", "cage");
						fprintf(F_out_type, "%s,", "cage");
						fprintf(f_out_dl, "%s,", "cage");
						fprintf(f_liste, "%s,", "cage");
					}
					else if(classification[i].classe == 2)
					{
						fprintf(F_out, "%s,", "pre cage");
						fprintf(F_out_type, "%s,", "pre cage");
						fprintf(f_out_dl, "%s,", "pre cage");
						fprintf(f_liste, "%s,", "pre cage");
					}
					else
					{
						fprintf(F_out, "%s,", "non cage");
						fprintf(F_out_type, "%s,", "non cage");
						fprintf(f_out_dl, "%s,", "non cage");
						fprintf(f_liste, "%s,", "non cage");
					}
				}
			}
			
			/*
			 *  Construction graphe cycle
			 * 
			 * */
			init_atom_num(); // associe les numéros atomiques aux symboles chimiques
			struct molecule m = lire_molecule_mol(F); // renvoie un graphe moléculaire à partir du fichier .mol
			//afficheInfosGrapheMol(m);
			GRAPHE_CYCLE cy = construction_graphe_cycles(m); // comme défini dans la thèse de Stefi (avec l'union des bases de cycles comme ensemble de sommets)
			cy = remove_edges(cy); // supprime les arêtes de type 3 (chaîne entre 2 cycles, voir Thèse Stefi) et les arêtes de poids == 0 (cad juste un sommet en commun)
			
			//printf("nb sommets: %d nb_aretes:%d\n", cy.nb_sommets, cy.nb_aretes);
		
			//afficheInfosGrapheCycle(cy);
			
			/*
			 * 
			 *  Calcul des cliques de taille 3 
			 * 
			 * */
			int *degre = calcul_degre(cy);
			int *store = malloc(cy.nb_sommets * sizeof(int));
			int n; // nb de bouboules (= coins artificiels)
			int c; // nb de coins qu'on garde
			int** tab_bouboule = init_tab_deux_dimensions(NB_TAB, NB_TAB);
			
			
			// d'abord on cherche les cliques pour détecter les bouboules
			n = findCliquesToDetectBouboules(cy, 0, 1, 3, degre, store, tab_bouboule, 0);
			// puis on recherche les cliques pour les compter par type 
			/* type 0: bouboule (coin artifciel)
			 * type 1: coin bouboule (coin quasi-artificiel)
			 * type 2: coin pas ouvert
			 * type 4: coin avec double liaisons aux jonctions
			 * type 3: tous les autres
			 * */
			c = findCliquesToCount(cy, 0, 1, 3, degre, store, F_out, F_out_type, f_out_dl, tab_bouboule, n, 0, m);
			fprintf(F_out, ",");
			
			/*
			 * 
			 *  Construction du graphe de coins
			 * 
			 * */
			GRAPHE_COIN gc;
			gc = construct_graph_coin(cy, c, degre, store, tab_bouboule, n, m);
			
			//afficheInfosGrapheCycle(cy);
			//afficheInfosGrapheCoin(gc);
			
			/*
			 * Génération des fichiers de visualisation
			 * 
			 * */
			
			genere_dot_file_cycles(cy, name, taille);
			genere_dot_file_coins(gc, name, taille);
			/*
			 * 
			 *  Lister les coins avec leurs caractéristiques (poids des sommets, poids des arêtes)
			 * */

			for(i=0;i<gc.nb_sommets;i++)
			{
				fprintf(f_liste, "%d%d%d%s%d%d%d,", gc.liste_sommets[i].poids[0], gc.liste_sommets[i].poids[1], gc.liste_sommets[i].poids[2], " ", gc.liste_sommets[i].liaisons_communs[0], gc.liste_sommets[i].liaisons_communs[1], gc.liste_sommets[i].liaisons_communs[2]);
			}
			fprintf(f_liste, "\n");
			
			// Désallocation
			free(store);
			free(degre);
			
			fprintf(F_out, "\n");
			fprintf(F_out_type,"\n");
			fprintf(f_out_dl, "\n");
			//printf("\n");
			
			fclose(F);
		
			strncpy(name, "", taille-4);
			free(name);
			free(name_file_mol);
			
			
			for(i=0;i<NB_TAB;i++)
			{
				free(tab_bouboule[i]);
			}
			free(tab_bouboule);
			
			if(taille_base > 0)
			{
				for( i = 0; i < taille_base;i++)
				{
					liberer_un_cycle(labase[i]);
				}
				free(labase);
				
			}
			
			liberer_graphe_cycles(cy);
			liberer_molecule(m);
			liberer_graphe_coin(gc);
			
			
		}
	}
	fclose(F_out);
	fclose(F_out_type);
	fclose(f_out_dl);
	fclose(f_liste);
	closedir(rep);
	free(classification);
	exit(0);
	
}


