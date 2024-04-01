#include "graphe_cycles.h"

void probleme_memoire()
{
	fprintf(stdout,"Cannot allocate memory\n");
		exit(145);
}

graphemol conversion_mol_graphe(struct molecule m)
{
	graphemol g;
	g.nb_atomes = m.nb_atomes;
	g.nb_liaisons = m.nb_liaisons;
	g.liste_atomes = malloc(g.nb_atomes *sizeof(int));
	g.type_atomes = malloc(g.nb_atomes *sizeof(int));
	g.liste_liaisons = malloc(g.nb_liaisons *sizeof(struct liaison));
	g.nb_connexe = 0;
	g.pere = 0;
	if(g.liste_atomes == NULL)
		probleme_memoire();
	if(g.type_atomes == NULL)
		probleme_memoire();
	
	int i,j;
	for (i = 0; i < g.nb_atomes; i++)
	{
		g.liste_atomes[i] = i + 1;
		g.type_atomes[i] = m.liste_atomes[i];
	}
	for (i = 0; i < g.nb_liaisons; i++)
	{
		g.liste_liaisons[i].A1 = m.liste_liaisons[i].A1;
		g.liste_liaisons[i].A2 = m.liste_liaisons[i].A2;
		g.liste_liaisons[i].l_type = m.liste_liaisons[i].l_type;
	}
	g.matrice_liaisons = malloc(g.nb_atomes *sizeof(int *));
	if( g.matrice_liaisons == NULL)
		probleme_memoire();
	for (i = 0; i < g.nb_atomes; i++)
	{
		g.matrice_liaisons[i]= malloc(g.nb_atomes *sizeof(int));
		if(g.matrice_liaisons[i] == NULL)
			probleme_memoire();
	}

	for (i = 0; i < g.nb_atomes; i++)
	{
		for (j = 0; j < g.nb_atomes; j++)
		{
			g.matrice_liaisons[i][j] = m.matrice_liaisons[i][j];
		}
	}
	//affiche_graphemol(g);
	return g;
}

int position_graphemol(graphemol g,int element)
{
	int i ;
	int pos = -1;
	for ( i = 0; i < g.nb_atomes;i++)
	{
		if(g.liste_atomes[i]== element)
		{	pos = i;
			break;
		}
	}
	if( pos == -1)
	{
		printf("%d\n", element);
		fprintf(stdout, "Impossible de trouver la position de cet atome dans la molecule\n");
		exit(458);
	}
	return pos;
}
int *calcul_degre_mol( graphemol m)

{
	int *degre = malloc( m.nb_atomes * sizeof( int)); // le degré de tous les sommets du graphe

	int i; 
	for(i = 0; i < m.nb_atomes; i++)
		degre[i] = 0;
	//remplissage
	int pos1,pos2;
	for ( i = 0; i < m.nb_liaisons; i++)
	{
		pos1 = position_graphemol(m,m.liste_liaisons[i].A1);
		pos2 = position_graphemol(m,m.liste_liaisons[i].A2);
		degre[pos1]++;
		degre[pos2]++;
	}

	return degre;
}


void liberer_graphemol(graphemol g)
{
	
	if ( g.liste_atomes != NULL)
		free(g.liste_atomes);
	if ( g.type_atomes != NULL)
		free(g.type_atomes);
	if( g.liste_liaisons != NULL)
		free(g.liste_liaisons);
	
	int i;
	for( i = 0; i < g.nb_atomes;i++)
	{
		if(g.matrice_liaisons[i] != NULL)
			free(g.matrice_liaisons[i]);
	}
	if(g.matrice_liaisons != NULL)
		free(g.matrice_liaisons);
	

}
graphemol suppression_aretes_mol(graphemol m, int i,int *degre)
{
	graphemol M;
	if(degre[i - 1] == 1)
		M.nb_liaisons = m.nb_liaisons - 1;
	else
		M.nb_liaisons = m.nb_liaisons;
	M.nb_atomes = m.nb_atomes;
	M.liste_atomes   = malloc(M.nb_atomes   * sizeof(int));
	M.type_atomes   = malloc(M.nb_atomes   * sizeof(int));
	M.liste_liaisons = malloc(M.nb_liaisons * sizeof(struct liaison));
	M.matrice_liaisons = malloc(M.nb_atomes   * sizeof(int *));
	M.nb_connexe = m.nb_connexe;
	M.pere = m.pere;
	int t;
	for ( t = 0; t < m.nb_atomes;t++)
		M.matrice_liaisons[t] = malloc(M.nb_atomes *sizeof(int));
	int a,k = 0; 
	//mise a jour des atomes
	for (a = 0; a < m.nb_atomes ; a++) 
	{
		
			M.type_atomes[k] = m.type_atomes[a];
			M.liste_atomes[k] = m.liste_atomes[a];
			k++;
	}	
	for ( t = 0; t < m.nb_atomes;t++)
	{
		for ( a = 0; a < M.nb_atomes; a++)
			M.matrice_liaisons[t][a] = 0;
	}
	k = 0;
	for (a = 0; a < m.nb_liaisons; a++) 
	{
		if(m.liste_liaisons[a].A1 != i && m.liste_liaisons[a].A2 != i )
		{
			M.liste_liaisons[k].A1 = m.liste_liaisons[a].A1;
			M.liste_liaisons[k].A2 = m.liste_liaisons[a].A2 ;
			M.liste_liaisons[k].l_type = m.liste_liaisons[a].l_type;
			M.matrice_liaisons[M.liste_liaisons[k].A1 - 1][M.liste_liaisons[k].A2 -1] = M.liste_liaisons[k].l_type;
			k++;
		}
	}	
	liberer_graphemol(m);
	return M;
}

// supprime les sommets de degré 0 ou 1 dans g
graphemol modification_structure_mol(graphemol g, int* degre,int* deja_elimine)
{
	int i = verif_mol(g,degre,deja_elimine) ;
	//affichage_degre(m,degre); 
	while(i != -1)
	{
		//fprintf(stdout, "le sommet %d\n",i +1  );
		g = suppression_aretes_mol( g, i + 1,degre);

		if ( degre != NULL)
			free(degre);
		degre = calcul_degre_mol(g);
		i = verif_mol(g,degre,deja_elimine) ; 
	}
	free(degre);
	return g;
}

// renvoie le premier sommet de degré 0 ou 1
int verif_mol(graphemol g, int* degre, int *deja_elimine)
{
	int i ,verif = -1; 

	for ( i = 0; i < g.nb_atomes; i++)
	{
		if (degre[i] == 1 || degre[i] == 0 )
		{
			if(deja_elimine[i] == 0)
			{
				verif = i;
				deja_elimine[i] = 1;
				break;
			}
			
		}
			
	}
	return verif;
}

struct liste_voisins*  construction_voisinage_graphemol(graphemol M){
	
	struct liste_voisins *voisins;
	int i,j,k,a1, a2;
	int nb_sommets = M.nb_atomes;
	int tableau_nb_voisins[nb_sommets];
	//initialisation du tableau à 0
	for(i = 0; i < nb_sommets; i++)
		tableau_nb_voisins[i] = 0;
	//remplissage
	int pos1,pos2;	
	for(i = 0; i < M.nb_liaisons; i++)
	{
		
		a1 = M.liste_liaisons[i].A1;
		a2 = M.liste_liaisons[i].A2;
		pos1= position_graphemol(M,a1);
		pos2= position_graphemol(M,a2);
		tableau_nb_voisins[ pos1]++;
		tableau_nb_voisins[ pos2]++;
		
	
	}
	//allocation memoire
	voisins = malloc( nb_sommets * sizeof(struct liste_voisins));
	if( voisins == NULL)
		probleme_memoire();
	
	for(i = 0; i < nb_sommets; i++)
	{
		//printf("i = %d; atome = %d et voisins = %d \n",i,M.liste_atomes[i], tableau_nb_voisins[i]);
		voisins[i].id_atome = M.liste_atomes[i] ;
		voisins[i].nb_voisins = tableau_nb_voisins[i];
		if( voisins[i].nb_voisins  > 0 )
		{
		
			voisins[i].id_voisins = malloc(voisins[i].nb_voisins * sizeof(int));
			
			if(voisins[i].id_voisins == NULL)
				probleme_memoire();
			k = 0;
		
			for(j = 0; j < M.nb_liaisons; j++)
			{
				if( M.liste_liaisons[j].A1 == voisins[i].id_atome )
				{
					voisins[i].id_voisins[k] = M.liste_liaisons[j].A2;
					//printf("%d %d\n", k+1,voisins[i].id_voisins[k]);
					k++;
				}
				
				else if ( M.liste_liaisons[j].A2 == voisins[i].id_atome )
				{
					voisins[i].id_voisins[k] = M.liste_liaisons[j].A1;
					//printf("%d %d\n", k+1,voisins[i].id_voisins[k]);
					k++;
				}
				
				
			}
		} 

	}
	//affichage_liste_voisinage_graphemol(voisins,M);
	return voisins;
}

//trouve tous les sommets atteignables a partir du sommet i

int * sommets_atteignable_graphemol(int i,graphemol m,struct liste_voisins * v,int *deja_elimine)
{
	int k,l,somme,somme_new,nb_valeurs,pos;
	int *trouve = malloc(m.nb_atomes* sizeof(int));

	if( trouve == NULL)
		probleme_memoire();
	
	//initialise à distance non atteignable pour tous les autres atomes
	for (k = 0; k < m.nb_atomes; k++)
		trouve[k] = -1;
	//metre à jour les voisins les plus proches ( distance 1)
	for (k = 0; k < v[i].nb_voisins; k++)
	{
		pos = position_graphemol(m,v[i].id_voisins[k]);
		trouve[pos] = 1;
	}
		
	trouve[i] = 0;
	//parcours en largeur
	int total = 0;
	for ( k = 0; k < m.nb_atomes;k++)
	{
		if(deja_elimine[i] == 1)
			total++;
	}
	somme = (v[i].nb_voisins - m.nb_atomes + total);  //nombre de d'atomes pas encore atteint
	somme_new = 0;
	int iteration = 1;
	while (somme != somme_new )
	{
		somme_new = somme;
		nb_valeurs = 0;
		for (k = 0; k < m.nb_atomes; k++)
		{
			if(trouve[k] == iteration)
			{
				//printf("iteration = %d sommet = %d\n",iteration ,k +1);
				for (l = 0; l < v[k].nb_voisins; l++)
				{
					pos = position_graphemol(m,v[k].id_voisins[l]);
					if( trouve[pos ] == -1)
					{
						trouve[pos] = trouve[k] + 1;
						nb_valeurs ++;
					}
						
				}
			}
		}
		somme = somme_new + nb_valeurs;
		iteration = iteration + 1;
	}
	return trouve;

}

//liberation liste_voisins


void liberation_liste_voisins_graphemol(struct liste_voisins *v,graphemol m)
{
	int i;
//liberation memoire liste des voisins 
	for (i = 0; i < m.nb_atomes; i++)
	{
		if(v[i].nb_voisins > 0)
			free(v[i].id_voisins);
	}
	free(v);
}

//construit tous les sous graphe induits connexes de la molecule
graphemol * ensemble_connexe_graphemol(graphemol m , int *deja_elimine)
{
	graphemol *g = NULL;
	int *degre = calcul_degre_mol(m);
	int marque[m.nb_atomes];
	int i,j,position,somme = m.nb_atomes;
	for ( i = 0; i < m.nb_atomes; i++)
	{
		if( deja_elimine[i] == 1)
			somme--;
	}

	//affiche_mol(m);
	int nb_connexe = 0;
	struct liste_voisins *v =  construction_voisinage_graphemol(m); // liste des voisins de chaque atome 
	for ( i = 0; i < m.nb_atomes; i++)
	{
		if(deja_elimine[i] == 0)
			marque[i] = 1;
		else
			marque[i] = 0;
	}
	
	// calcul des composantes connexes du graphe moléculaire (sans les atomes qui n'appartiennent pas à des cycles)
	int *t;
	int sommets;
	while( somme > 0)
	{
		nb_connexe++;
		for ( i = 0; i < m.nb_atomes; i++)
		{
			if(marque[i] == 1)
			{
				position = i;
				break;
			}
		}
		//printf("position = %d \n", position);
		t  = sommets_atteignable_graphemol(position,m,v,deja_elimine); // ensemble des sommets atteignables à partit de l'atome position (et la distance)

		sommets = 0;
		for ( i = 0; i < m.nb_atomes; i++)
		{	
			//printf("%d ",t[i] );
			if(t[i] >= 0 && marque[i] == 1)
			{
				marque[i] = 0;
				sommets = sommets + 1;
			}

		}
		
		if (g == NULL && nb_connexe == 1)
		{
			//printf("nb connexes %d\n",nb_connexe );
			g = malloc(nb_connexe * sizeof(graphemol));
			if( g == NULL)
				probleme_memoire();


		
		}
		else
		{	
			//printf("nb connexes %d\n",nb_connexe );
			g = realloc(g,nb_connexe * sizeof(graphemol));
			if( g == NULL)
				probleme_memoire();
		}
		
		//mise a jour de la composante 
		g[nb_connexe - 1].nb_atomes = sommets;
		g[nb_connexe - 1].nb_liaisons = 0;
		g[nb_connexe - 1].pere = m.pere + 1;
		int aretes = 0;

		//compter les aretes
		int pos1,pos2;
		for ( i = 0; i < m.nb_liaisons;i++)
		{
			pos1 = position_graphemol(m,m.liste_liaisons[i].A1);
			pos2 = position_graphemol(m,m.liste_liaisons[i].A2);
			if(t[pos1] >= 0 || t[pos2] >= 0)
				aretes++;
		}
		//stockage des sommets
		//printf("connexe , %d\n",nb_connexe );
		g[nb_connexe - 1].liste_atomes = malloc(sommets *sizeof(int));
		g[nb_connexe - 1].type_atomes = malloc(sommets *sizeof(int));

		if(g[nb_connexe - 1].liste_atomes== NULL)
			probleme_memoire();
		if(g[nb_connexe - 1].type_atomes== NULL)
			probleme_memoire();
		i = 0;
		int compteur = 0;
		while(i < m.nb_atomes)
		{
			
			if(t[i] >=0)
			{
				
				g[nb_connexe - 1].liste_atomes[compteur] = m.liste_atomes[i];
				g[nb_connexe - 1].type_atomes[compteur] = m.type_atomes[i];
				compteur++;
			}
			
			i = i + 1 ;
		}
		//stockage des aretes
		g[nb_connexe - 1].nb_liaisons = aretes; 
		g[nb_connexe - 1].liste_liaisons = malloc(aretes *sizeof(struct liaison));
		if(g[nb_connexe - 1].liste_liaisons == NULL)
			probleme_memoire();
		//stockage du type des liaisons
		g[nb_connexe - 1].nb_connexe = nb_connexe;
		g[nb_connexe - 1].matrice_liaisons = malloc(g[nb_connexe - 1].nb_atomes *sizeof(int *));
		for ( i = 0; i < g[nb_connexe - 1].nb_atomes; i++)
		{
			g[nb_connexe - 1].matrice_liaisons[i] = malloc(g[nb_connexe - 1].nb_atomes * sizeof(int));
			if(g[nb_connexe - 1].matrice_liaisons[i] == NULL)
				probleme_memoire();	
		}
		

		for ( i = 0; i < g[nb_connexe - 1].nb_atomes; i++)
		{
			for ( j = 0; j < g[nb_connexe - 1].nb_atomes; j++)
				g[nb_connexe - 1].matrice_liaisons[i][j] = 0;
		}


		i = 0;
		compteur = 0;
		while ( i < m.nb_liaisons)
		{
			pos1 = position_graphemol(m,m.liste_liaisons[i].A1);
			pos2 = position_graphemol(m,m.liste_liaisons[i].A2);
			if(t[pos1] >= 0 || t[pos2] >= 0)
			{
				g[nb_connexe - 1].liste_liaisons[compteur].A1 = m.liste_liaisons[i].A1;
				g[nb_connexe - 1].liste_liaisons[compteur].A2 = m.liste_liaisons[i].A2;
				g[nb_connexe - 1].liste_liaisons[compteur].l_type = m.liste_liaisons[i].l_type;
				compteur++;

			}
			i = i + 1 ;
		}
		for ( i = 0; i < g[nb_connexe -1].nb_liaisons; i++)
		{
			pos1 = position_graphemol( g[nb_connexe - 1],g[nb_connexe - 1].liste_liaisons[i].A1);
			pos2 = position_graphemol( g[nb_connexe - 1],g[nb_connexe - 1].liste_liaisons[i].A2);

			g[nb_connexe - 1].matrice_liaisons[pos1][pos2]=g[nb_connexe - 1].liste_liaisons[i].l_type;
			g[nb_connexe - 1].matrice_liaisons[pos2][pos1]=g[nb_connexe - 1].liste_liaisons[i].l_type;
		}

		



		//recalcul de la somme
		somme = 0;
		for ( i = 0; i < m.nb_atomes; i++)
			somme += marque[i];
		//printf("somme = %d\n", somme);

		free(t);
		
	}
	free(degre);
	liberation_liste_voisins_graphemol(v,m);

	g[0].nb_connexe = nb_connexe ;
	return g;
}

int nombre_isthmes_graphemol(graphemol g, int *deja_elimine)
{
	int nb_isthmes = 0,i;
	i = 0; 

	while ( i  < g.nb_liaisons)
	{
		i++;
		//printf("%d %d %d \n",g.liste_liaisons[0].A1,g.liste_liaisons[0].A2,nb_isthmes);
		nb_isthmes += retirer_liaison_connexe(g,g.liste_liaisons[0],deja_elimine);
	}


	return nb_isthmes;
}

graphemol enlever_une_arete(graphemol g, struct liaison l)
{
	int i,pos1,pos2;
	pos1 = position_graphemol(g,l.A1);
	pos2 = position_graphemol(g,l.A2);

	g.matrice_liaisons[pos1][pos2]= 0;
	g.matrice_liaisons[pos2][pos1]= 0;

	//stockage
	struct liaison store[g.nb_liaisons];
	for ( i = 0; i < g.nb_liaisons; i++)
	{
		store[i].A1 = g.liste_liaisons[i].A1;
		store[i].A2 = g.liste_liaisons[i].A2;
		store[i].l_type = g.liste_liaisons[i].l_type;
	}

	g.nb_liaisons = g.nb_liaisons - 1;

	int compteur = 0; 
	for ( i = 0; i < g.nb_liaisons + 1; i++)
	{	
		if(store[i].A1 != l.A1 || store[i].A2 != l.A2 )
		{
			g.liste_liaisons[compteur].A1 = store[i].A1;
			g.liste_liaisons[compteur].A2 = store[i].A2;
			g.liste_liaisons[compteur].l_type = store[i].l_type;
			compteur++;
		}
	}
	g.liste_liaisons[g.nb_liaisons ].A1 =0;
	g.liste_liaisons[g.nb_liaisons ].A2 =0;
	g.liste_liaisons[g.nb_liaisons ].l_type =0;

	return g;
}

//verifie si dans le graphe de la molecule m il existe une chaine entre les sommets en position i et j.( 1 si oui et 0 sinon)

int existe_chaine_graphemol(int i, int j , graphemol m,struct liste_voisins * v)
{
	int resultat = 0,k,l,somme,somme_new,nb_valeurs;
	int trouve[m.nb_atomes];
	int pos;
	//initialise à distance non atteignable pour tous les autres atomes
	for (k = 0; k < m.nb_atomes; k++)
		trouve[k] = -1;
	//metre à jour les voisins les plus proches ( distance 1)
	for (k = 0; k < v[i].nb_voisins; k++)
	{
		pos = position_graphemol(m,v[i].id_voisins[k]);
		trouve[pos] = 1;
	}
		
	trouve[i] = 0;
	//parcours en largeur
	somme = (v[i].nb_voisins - m.nb_atomes);  //nombre de d'atomes pas encore atteint
	somme_new = 0;
	int iteration = 1;
	while ( trouve[j] == -1 && somme != somme_new )
	{
		somme_new = somme;
		nb_valeurs = 0;
		for (k = 0; k < m.nb_atomes; k++)
		{
			if(trouve[k] == iteration)
			{
				//printf("iteration = %d sommet = %d\n",iteration ,k +1);
				for (l = 0; l < v[k].nb_voisins; l++)
				{
					pos = position_graphemol(m,v[k].id_voisins[l]);
					if( trouve[pos] == -1)
					{
						trouve[pos ] = trouve[k] + 1;
						nb_valeurs++;
					}
						
				}
			}
		}
		somme = somme_new + nb_valeurs;
		iteration = iteration + 1;
	}
	if( trouve[j] != -1)
		resultat = 1;

	return resultat;
}
int est_connexe_graphemol(graphemol m, int *deja_elimine)
{
	//affiche_mol(m);
	int resultat  = 1,i,j,chaine;
	//construction de la liste des voisins
	struct liste_voisins *v =  construction_voisinage_graphemol(m);

	for (i = 0; i < m.nb_atomes - 1; i++)
	{
		//trouver s'il existe un chemin de l'atome en position i vers l'atome en position jj
		for (j = i + 1; j < m.nb_atomes; j++)
		{
			//verifie s'il existe une chaine entre i et j dans m
			//printf(" c %d %d %d %d \n", m.liste_atomes[i],m.liste_atomes[j],deja_elimine[i],deja_elimine[j]);
			if(deja_elimine[i] == 0 && deja_elimine[j] == 0)
			{
				chaine = existe_chaine_graphemol(i,j,m,v);
				//printf("existe chaine %d %d %d\n", m.liste_atomes[i],m.liste_atomes[j],chaine);
				if(!chaine)
				{
					resultat = 0;
					break;
				}

			}
			
		}
		if(j != m.nb_atomes && !chaine)
			break;
	}
	for (i = 0; i < m.nb_atomes; i++)
	{
		if(v[i].nb_voisins > 0)
			free(v[i].id_voisins);
	}
	free(v);

	return resultat;
}


graphemol ajouter_une_arete(graphemol g, struct liaison l)
{
	int i,p,pos1=-1,pos2=-1;
	for( p = 0; p < g.nb_atomes;p++)
	{
		if(g.liste_atomes[p] == l.A1)
		{
			pos1=p;
			break;
		}
	}

	for( p = 0;p < g.nb_atomes;p++)
	{
		if(g.liste_atomes[p] == l.A2)
		{
			pos2=p;
			break;
		}
	}
	if(pos1 == -1 || pos2 ==-1)
	{
		fprintf(stdout, "La position de l'atome dans la molécule n'a pas été trouvé\n");
		exit(234);
	}	

	g.matrice_liaisons[pos1][pos2]= l.l_type;
	g.matrice_liaisons[pos2][pos1]= l.l_type;
	//stockage
	g.nb_liaisons = g.nb_liaisons + 1;
	for ( i = 0; i < g.nb_liaisons ; i++)
	{
		//printf("--%d %d %d ",i, g.liste_liaisons[i].A1,g.liste_liaisons[i].A2);
		if(g.liste_liaisons[i].A1 == 0 && g.liste_liaisons[i].A2 == 0)
		{
			g.liste_liaisons[i ].A1 =l.A1;
			g.liste_liaisons[i].A2 =l.A2;
			g.liste_liaisons[i ].l_type =l.l_type;
		}
	}
	
	

	return g;
}
int retirer_liaison_connexe(graphemol g, struct liaison l, int *deja_elimine)//on retire la liaison l et on verifie que le graphe reste connexe
{//on retourne 0 si le graphe reste connexe et 1 sinon
	
	int resultat = 0;


	g = enlever_une_arete(g,l);
	
	resultat = 1 - est_connexe_graphemol(g,deja_elimine);
	g = ajouter_une_arete(g,l);
	
	

	return resultat;
}

isthmes *retrouver_tous_isthmes(graphemol *liste_connexe,isthmes *liste, int *deja_elimine)
{
	int nb_connexe = liste_connexe[0].nb_connexe;
	isthmes p;
	int nb_isthmes = 0,i,j,l;
	int *elimine;
	for ( j = 0;  j < nb_connexe ;  j++)
	{
		i = 0; 
		//printf("j = %d nb liaisons = %d\n",j, liste_connexe[j].nb_liaisons );

		elimine = malloc(liste_connexe[j].nb_atomes * sizeof(int));
		for ( l = 0;  l < liste_connexe[j].nb_atomes;l++)
		{	
				elimine[l] = 0;
		}
		while ( i  < liste_connexe[j].nb_liaisons)
		{
			p.l.A1 = liste_connexe[j].liste_liaisons[0].A1;
			p.l.A2 = liste_connexe[j].liste_liaisons[0].A2;
			p.l.l_type = liste_connexe[j].liste_liaisons[0].l_type;
			p.id_composant = liste_connexe[j].pere;
			//printf("j: %d %d %d\n",liste_connexe[j].liste_liaisons[0].A1,liste_connexe[j].liste_liaisons[0].A2, retirer_liaison_connexe(liste_connexe[j],liste_connexe[j].liste_liaisons[0],deja_elimine));
			if(retirer_liaison_connexe(liste_connexe[j],liste_connexe[j].liste_liaisons[0],elimine) == 1)
			{
				liste[nb_isthmes].l.A1 = p.l.A1;
				liste[nb_isthmes].l.A2 = p.l.A2;
				liste[nb_isthmes].l.l_type = p.l.l_type;
				liste[nb_isthmes].id_composant = p.id_composant;
				//printf("i: %d %d \n",liste[nb_isthmes].l.A1,liste[nb_isthmes].l.A2 );  
				nb_isthmes++;
			}
			i++;
		}
		free(elimine);
	}
	//affiche_liste_isthmes(liste,nb_isthmes);
	return liste;
}


graphemol enlever_tous_isthmes(graphemol g,int *deja_elimine,isthmes *liste_isthmes,int nb_isthmes)
{

	graphemol h;
	int *degre = calcul_degre_mol(g);
	int atomes = 0,i,j,pos,pos1;
	for (i = 0; i < nb_isthmes; i++)
	{
		//printf("%d %d ---- %d %d\n", position_graphemol(g,liste_isthmes[i].l.A1),degre[position_graphemol(g,liste_isthmes[i].l.A1)],position_graphemol(g,liste_isthmes[i].l.A2),degre[position_graphemol(g,liste_isthmes[i].l.A2)]);
		pos = position_graphemol(g,liste_isthmes[i].l.A1);
		degre[pos]--;
		pos = position_graphemol(g,liste_isthmes[i].l.A2);
		degre[pos]--;
	}
	for(i = 0; i < g.nb_atomes; i++)
	{
		if(degre[i] > 0)
			atomes++;
	}
	h.nb_atomes = atomes;
	h.nb_liaisons = g.nb_liaisons - nb_isthmes;
	h.liste_atomes = malloc(h.nb_atomes * sizeof(int));
	if(h.liste_atomes == NULL)
		probleme_memoire();
	h.type_atomes = malloc(h.nb_atomes * sizeof(int));
	if(h.type_atomes == NULL)
		probleme_memoire();
	i = 0;

	while ( i < h.nb_atomes)
	{
		for(j = 0; j <g.nb_atomes;j++ )
		{
			if(degre[j] > 0)
			{
				h.liste_atomes[i] = g.liste_atomes[j];
				h.type_atomes[i] = g.type_atomes[j];
				i++;
			}
		}
	}
	h.liste_liaisons = malloc(h.nb_liaisons *sizeof(struct liaison));
	int compteur = 0 ;
	for( i  = 0; i < nb_isthmes; i++)
	{
		for(j = 0;  j < g.nb_liaisons;j++)
		{
			if((g.liste_liaisons[j].A1 == liste_isthmes[i].l.A1 )&&(g.liste_liaisons[j].A2 == liste_isthmes[i].l.A2))
			{
				g.liste_liaisons[j].A1 = 0;
				g.liste_liaisons[j].A2 = 0;
				break;
			}
		}
	}
	for(j = 0;  j < g.nb_liaisons;j++)
	{
		if((g.liste_liaisons[j].A1 != 0 )&&(g.liste_liaisons[j].A2 != 0))
		{
			h.liste_liaisons[compteur].A1 = g.liste_liaisons[j].A1 ;
			h.liste_liaisons[compteur].A2 = g.liste_liaisons[j].A2;
			h.liste_liaisons[compteur].l_type = g.liste_liaisons[j].l_type;
			compteur++;
		}
	}
	h.matrice_liaisons = malloc(h.nb_atomes * sizeof(int *));
	
	for ( i = 0; i < h.nb_atomes; i++)
		h.matrice_liaisons[i] = malloc(h.nb_atomes  *sizeof(int));
	for ( i = 0; i < h.nb_atomes; i++)
	{
		for ( j = 0; j < h.nb_atomes; j++)
		{
			h.matrice_liaisons[i][j] = 0;
		}
	}

	for(j = 0;  j < h.nb_liaisons;j++)
	{
		pos = position_graphemol(h,h.liste_liaisons[j].A1);
		pos1 = position_graphemol(h,h.liste_liaisons[j].A2);
		h.matrice_liaisons[pos][pos1] = h.liste_liaisons[j].l_type;
		h.matrice_liaisons[pos1][pos] = h.liste_liaisons[j].l_type;
	}
	h.nb_connexe  = g.nb_connexe;
	h.pere  = g.pere + 1;
	free(degre);
	return h;
}

int sommet_dans_basniveau(graphemol t,int p)
{
	int res = 0;
	int i;
	for(i = 0; i < t.nb_atomes;i++)
	{
		if(t.liste_atomes[i] == p)
		{
			res  = 1;
			break;
		}
	}

	return res;
}
//liberer la memoire de liste voisins 

void liberer_memoire_voisins(struct liste_voisins *v,graphemol m)
{
	int i;
	for (i = 0; i < m.nb_atomes; i++)
	{
		if(v[i].nb_voisins > 0)
			free(v[i].id_voisins);
	}
	free(v);

}
//calcule un plus court chemin entre deux sommets 

int *plus_court_chemin(int sommet1,int sommet2,graphemol m)
{

	int *chemin =  NULL;
	struct liste_voisins *v = construction_voisinage_graphemol(m);

	//affichage_liste_voisinage_graphemol(voisins,m);
	int distance[m.nb_atomes];
	int predecesseur[m.nb_atomes];

	int i,j;
	
	for (i = 0; i < m.nb_atomes; i++)
	{
		distance[i] = -1;
		predecesseur[i] = -1;
	}
		
	distance[position_graphemol(m,sommet1)] = 0;
	predecesseur[position_graphemol(m,sommet1)] = position_graphemol(m,sommet1);

	int pos = position_graphemol(m,sommet2);
	int iteration = 0;
	while(distance[pos] == -1)
	{
		for( i = 0; i < m.nb_atomes; i++)
		{
			if(distance[i] == iteration)
			{
				for(j = 0; j < v[i].nb_voisins;j++)
				{
					if(distance[position_graphemol(m,v[i].id_voisins[j])] == -1)
					{
						distance[position_graphemol(m,v[i].id_voisins[j])] = iteration + 1;
						predecesseur[position_graphemol(m,v[i].id_voisins[j])] = m.liste_atomes[i];
					}
				}
			}
		}
		iteration++;
	}
	chemin = malloc((distance[pos] +1) * sizeof(int));
	int s  = 0;
	chemin[s] = distance[pos];
	int pred = predecesseur[pos];
	s++;
	for(i = 0; i < distance[pos] - 1;i++)
	{
		chemin[s] = pred;
		pred  = predecesseur[position_graphemol(m,pred)];
		s++;
	}
	liberer_memoire_voisins(v,m);
	return chemin;
}

//calcule un plus court chemin entre deux sommets (avec que des sommets précédant sommet1 dans l'ordre des sommets)
int *plus_court_chemin_precedant(int sommet1,int sommet2,graphemol m)
{
	int *chemin =  NULL;

	// deja, si sommet 2 est après sommet 1 dans l'ordre qu'on a défini, on peut sortir
	if(position_graphemol(m, sommet1) <= position_graphemol(m, sommet2)) return NULL;
	

	struct liste_voisins *v = construction_voisinage_graphemol(m);
	//affichage_liste_voisinage_graphemol(voisins,m);
	int distance[m.nb_atomes];
	int predecesseur[m.nb_atomes];

	int mini, temp, pos_min;
	int i,j,k;
	
	
	// ordonner le voisinage
	for(i=0;i<m.nb_atomes;i++)
	{
		for(k=0;k<v[i].nb_voisins;k++)
		{
			mini = m.nb_atomes+1;
			pos_min = k;
			for(j=k;j<v[i].nb_voisins;j++)
			{
				if(position_graphemol(m, v[i].id_voisins[j]) < mini) 
				{
					mini = position_graphemol(m, v[i].id_voisins[j]);
					pos_min = j;
				}
			}
			if(pos_min != k)
			{
				temp = v[i].id_voisins[pos_min];
				v[i].id_voisins[pos_min] = v[i].id_voisins[k];
				v[i].id_voisins[k] = temp;
			}
		}
	}
	
	for (i = 0; i < m.nb_atomes; i++)
	{
		distance[i] = -1;
		predecesseur[i] = -1;
	}
		
	distance[position_graphemol(m,sommet1)] = 0;
	predecesseur[position_graphemol(m,sommet1)] = position_graphemol(m,sommet1);

	int pos = position_graphemol(m,sommet2);
	int iteration = 0;
	while(distance[pos] == -1)
	{
		for( i = 0; i < m.nb_atomes; i++)
		{
			if(distance[i] == iteration)
			{
				for(j = 0; j < v[i].nb_voisins;j++)
				{
					if(position_graphemol(m, v[i].id_voisins[j]) < position_graphemol(m, sommet1) && distance[position_graphemol(m,v[i].id_voisins[j])] == -1) // on ne modifie la distance que si le voisin est avant le sommet 1 dans l'ordre des sommets
					{
						//printf("i: %d\n",m.liste_atomes[i]);
						distance[position_graphemol(m,v[i].id_voisins[j])] = iteration + 1;
						predecesseur[position_graphemol(m,v[i].id_voisins[j])] = m.liste_atomes[i];
					}
				}
			}
		}
		iteration++;
	}

	
	chemin = malloc((distance[pos] +1) * sizeof(int));
	int s  = 0;
	chemin[s] = distance[pos];
	int pred = predecesseur[pos];
	s++;

	for(i = 0; i < distance[pos] - 1;i++)
	{
		if(position_graphemol(m, pred) >= position_graphemol(m, sommet1))
		{
			liberer_memoire_voisins(v,m);
			free(chemin);
			return NULL;
		}
		chemin[s] = pred;
		pred  = predecesseur[position_graphemol(m,pred)];
		s++;
	}

	
	liberer_memoire_voisins(v,m);
	return chemin;
}


// retourne 1 si les chemins sont independants et 0 sinon 
int chemin_independant(int *chemin1, int *chemin2,struct liaison l)
{
	
	int resultat = 1,trouve = 0;
	int i,j;
	for( i = 0; i < chemin1[0] -1;i++)
	{
		if(chemin1[i+1] == l.A1 ||chemin1[i+1] == l.A2)
		{
			resultat = 0;
			return resultat;
		}
	}

	for( i = 0; i < chemin2[0] -1;i++)
	{
		if(chemin2[i+1] == l.A1 || chemin2[i+1] == l.A2)
		{
			resultat = 0;
			return resultat;
		}
	}
	if(chemin1[0] > chemin2[0])
	{

		for( i = 0 ; i < chemin2[0] - 1; i++)
		{
			for( j = 0; j < chemin1[0] - 1; j++)
			{
				if(chemin1[j + 1] == chemin2[i + 1])
				{
					resultat = 0;
					trouve = 1;
					break;
				}
			}
			if( trouve  == 1)
				break;
		}
	}
	else
	{
		for( i = 0 ; i < chemin1[0] - 1; i++)
		{
			for( j = 0; j < chemin2[0] - 1; j++)
			{
				if(chemin1[i + 1] == chemin2[j + 1])
				{
					resultat = 0;
					trouve = 1;
					break;
				}
			}
			if( trouve  == 1)
				break;
		}

	}
	return resultat;
}

//cree le cycle qui passe le sommet et l'arete dans le graphe m
cycles creer_un_cycle(graphemol m , int sommet, struct liaison l, int *chemin1,int *chemin2)
{
	cycles c;
	int tableau[m.nb_atomes],i;

	for (i = 0; i < m.nb_atomes; i++)
		tableau[i] = 0;

	tableau[position_graphemol(m,sommet)] = 1;
	tableau[position_graphemol(m,l.A1)] = 1;
	tableau[position_graphemol(m,l.A2)] = 1;
	for(i  = 0; i < chemin1[0] - 1; i++)
		tableau[position_graphemol(m, chemin1[i+1])] = 1;
	for(i  = 0; i < chemin2[0] - 1; i++)
		tableau[position_graphemol(m, chemin2[i+1])] = 1;

	c.nb_atomes = 0;

	for (i = 0; i < m.nb_atomes; i++)
		c.nb_atomes+= tableau[i] ;

	c.sommet = sommet;
	c.arete  = l;
	c.chemin1 = malloc(chemin1[0] * sizeof(int));
	for( i = 0; i < chemin1[0]; i++)
		c.chemin1[i] = chemin1[i];
	
	c.chemin2 = malloc(chemin2[0] * sizeof(int));
	for( i = 0; i < chemin2[0]; i++)
		c.chemin2[i] = chemin2[i];
	
	c.sommets = malloc(c.nb_atomes * sizeof(int));
	int k = 0;
	for ( i = 0; i < m.nb_atomes; i++)
	{
		if(tableau[i] == 1)
		{
			c.sommets[k] = m.liste_atomes[i];
			k++;
		}
	}
	c.pere = 0;
	c.id_cycle = 0;
	return c;
}

//verifie si le cycle n'a pas encore été trouvé. retourne 1 si le cycle est deja das la base et 0 sinon
int verification_ajout_cycle(cycles *l, int nb_cycles , cycles c,graphemol m )
{
	int trouve = 0,i,j;
	
	int compteur ;	
	for( i =  0; i < nb_cycles; i++)
	{
		
		if(l[i].nb_atomes == c.nb_atomes)
		{
			compteur = 0;
			for( j =  0; j < c.nb_atomes; j++)
			{
				if(c.sommets[j] == l[i].sommets[j])
				{
					compteur++;
				}
			}
			if( compteur == c.nb_atomes)
			{
				trouve  = 1;
				break;
			}
		}
	}
	return trouve;
}

cycles *ajouter_un_cycle(cycles *liste, int nb_cycles , cycles c,graphemol m)
{
	//afficher_un_cycle(c,m);
	if( liste == NULL)
	{
		liste = malloc((nb_cycles + 1)* sizeof(cycles));

	}
	else
	{
		liste = realloc(liste, (nb_cycles + 1)* sizeof(cycles));
	}
	if( liste == NULL)
			probleme_memoire();
	liste[nb_cycles].nb_atomes = c.nb_atomes;
	liste[nb_cycles].sommet = c.sommet;
	liste[nb_cycles].pere = c.pere;
	liste[nb_cycles].arete.A1 = c.arete.A1;
	liste[nb_cycles].arete.A2 = c.arete.A2;
	liste[nb_cycles].sommets = malloc(c.nb_atomes * sizeof(int));
	if(liste[nb_cycles].sommets ==  NULL)
		probleme_memoire();
	int i;
	for( i = 0; i < c.nb_atomes; i++)
		liste[nb_cycles].sommets[i] = c.sommets[i];
	liste[nb_cycles].chemin1 = malloc(c.chemin1[0] *sizeof(int));
	if(liste[nb_cycles].chemin1 == NULL)
		probleme_memoire();
	for( i = 0; i < c.chemin1[0]; i++)
		liste[nb_cycles].chemin1[i] = c.chemin1[i];

	liste[nb_cycles].chemin2 = malloc(c.chemin2[0]  *sizeof(int));
	if(liste[nb_cycles].chemin2 == NULL)
		probleme_memoire();
	for( i = 0; i < c.chemin2[0] ; i++)
		liste[nb_cycles].chemin2[i] = c.chemin2[i];
	liste[nb_cycles].id_cycle = nb_cycles + 1;
	return liste;
}


//liberer un cycle

void liberer_un_cycle(cycles c)
{

	if( c.sommets != NULL)
		free(c.sommets);
	if(c.chemin1 != NULL)
		free(c.chemin1);
	if(c.chemin2 != NULL)
		free(c.chemin2);
}

int *concatener_deux_chemins(cycles c)
{
	int s = c.chemin1[0] +  c.chemin2[0] + 1;


	int *chemin  = malloc((s+1) * sizeof(int));
	chemin[0] = c.sommet;
	int i,pos = 1;
	for(i = c.chemin1[0] - 1; i > 0; i--)
	{
		chemin[pos] = c.chemin1[i];
		//printf("%d ", chemin[pos]);
		pos++;
	}
	chemin[pos] = c.arete.A1;
	
	pos++;
	chemin[pos] = c.arete.A2;
	pos++;
	for(i = 0; i < c.chemin2[0] - 1; i++)
	{
		
		chemin[pos] = c.chemin2[i + 1];
		pos++;
	}

	chemin[s] = c.sommet;
	return chemin;
}

int position_de_arete(int sommet1, int sommet2,graphemol m)
{
	int pos = -1;
	int i;

	for ( i = 0; i < m.nb_liaisons; i++)
	{
		if((m.liste_liaisons[i].A1 == sommet1 && m.liste_liaisons[i].A2 == sommet2 )||(m.liste_liaisons[i].A2 == sommet1 && m.liste_liaisons[i].A1 == sommet2 ))
		{
			pos = i;
			break;
		}
	}

	if( pos == -1)
	{
		fprintf(stdout, "impossible de trouver l'arete [ %d - %d] dans le graphe\n", sommet1,sommet2);
		exit(600);
	}
	return pos;
}


int fonction_xor(int a , int b)
{
	int c;
	if( a == 0)
	{
		if( b ==0)
			c = 0;
		else
			c = 1;

	}
	else
	{
		if( b == 0)
			c = 1;
		else 
			c = 0;
	}
	return c;
}


int *produit_xor_matrice(int *t, int *tab1, int *tab2,int nb)
{
	t= malloc(nb*sizeof(int));
	int i;
	for ( i = 0; i < nb; i++)
		t[i] = fonction_xor(tab1[i],tab2[i]);
	return t;
	
}
//verifie si deux tableaux sont les memes

int verification_egalite_tableaux(int *tab1,int *tab2 , int taille)
{
	
	int i, verif = 1 ;
	for( i = 0; i < taille ; i++)
	{
		if( tab1[i] != tab2[i])
		{
			verif = 0;
			break;
		}
	}
	
	return verif;
	
}

void ajouter_cycle_base(cycles c,graphemol m)

{
	//if (c.nb_atomes > taille_cycle)
	//	return; uniquement quand on limite la taille des cycles
	if(taille_base == 0)
	{
		labase = malloc((taille_base +1)* sizeof(cycles));
	}
	else
	{
		labase = realloc(labase,(taille_base +1)* sizeof(cycles));
	}
	if( labase == NULL)
			probleme_memoire();
	labase[taille_base] = creer_un_cycle(m,c.sommet,c.arete,c.chemin1,c.chemin2);
	labase[taille_base].id_cycle = taille_base;
	labase[taille_base].pere = c.pere;
	taille_base++;


}

int* init_t_zeros(int nb_liaisons)
{
	int *t = malloc(nb_liaisons*sizeof(int));
	int i;
	for(i=0;i<nb_liaisons;i++)
	{
		t[i] = 0;
	}
	return t;
}


int nb_1_in_matrix(int** matrice, int i, int nb_col)
{
	int j;
	int nb_1 = 0;
	for(j=0;j<nb_col;j++)
	{
		if(matrice[i][j] == 1) nb_1++;
		
	}
	return nb_1;
}


int* copie_chemin(int* chemin_vieux)
{
	//printf("taille chemin_vieux : %d\n", chemin_vieux[0]);
	int* chemin_new = malloc((chemin_vieux[0]+1)*sizeof(int));
	int i;
	
	for(i=0;i<chemin_vieux[0]+1;i++)
	{
		chemin_new[i] = chemin_vieux[i];
	}
	return chemin_new;
}

// génère la famille d'un cycle prototype (fonction list_paths article)
// on cherche tous les plus courts chemins de x à r en passant par des sommets de numéro plus petits que r et que x
int** list_path(int x, int r, int* chemin, struct graphe_dr g_dr, int* nb_chemins_courants, int taille_chemin, int** chemins_courants, int l_A)
{
	// on considère que dans chemin, on a le chemin à l'envers avec le début et la fin
	int i,j,k, nb_temp;
	int** chemins_courants_temp;
	int* chemin_t;
	if(x == r)
	{
		
		int** chemins_courants_t = malloc(1*sizeof(int*));
		chemins_courants_t[0] = copie_chemin(chemin);
		(*nb_chemins_courants) ++;

		return chemins_courants_t;
	}
	
	if(chemin == NULL)
	{
		chemin = malloc((taille_chemin)*sizeof(int));
		chemin[0] = 1;
		for(i=1;i<taille_chemin;i++)
		{
			chemin[i] = -1;
		}
	}
	else
	{
		chemin[0]++;
	}
	chemin[chemin[0]] = x;
	

	for(i=0;i<g_dr.nb_arcs;i++)
	{
		if(g_dr.liste_arcs[i].id1 == x)
		{
			nb_temp = (*nb_chemins_courants);
			chemin_t = copie_chemin(chemin);
		
			chemins_courants_temp = list_path(g_dr.liste_arcs[i].id2, r, chemin, g_dr, nb_chemins_courants, taille_chemin, chemins_courants, l_A);
			
			if(g_dr.liste_arcs[i].id2 == r) 
			{
				if(chemins_courants == NULL) chemins_courants = malloc((*nb_chemins_courants)*sizeof(int*));
				else chemins_courants = realloc(chemins_courants, (*nb_chemins_courants)*sizeof(int*));
				
				for(j=nb_temp;j<(*nb_chemins_courants);j++)
				{
					chemins_courants[j] = copie_chemin(chemins_courants_temp[j-nb_temp]);
					free(chemins_courants_temp[j-nb_temp]);
				}
				free(chemins_courants_temp);
				
			}
			else
			{
				chemins_courants = chemins_courants_temp;
			}
			for(j=0;j<chemin[0]+1;j++)
			{
				chemin[j] = chemin_t[j];
			}
			free(chemin_t);
		}
	}
	if(x == l_A) free(chemin); // dernier appel de list path, l_A est notre sommet de départ
	return chemins_courants;
}

void affichage_cycle(cycles c)
{
	printf("cycle : sommet %d, arete %d %d\n", c.sommet, c.arete.A1, c.arete.A2);
	int i;
	printf("chemin 1 : (taille %d) ",c.chemin1[0]);
	for(i=1;i<c.chemin1[0];i++)
	{
		printf("%d ", c.chemin1[i]);
	}
	printf("\n");
	printf("chemin 2 : (taille %d) ",c.chemin2[0]);
	for(i=1;i<c.chemin2[0];i++)
	{
		printf("%d ", c.chemin2[i]);
	}
	printf("\n");
}

int compte_1_ligne(int** matrice, int nb_colonnes, int num_ligne)
{
	int i;
	int compte = 0;
	for(i=0;i<nb_colonnes;i++)
	{
		if(matrice[num_ligne][i] == 1) compte ++;
	}
	return compte;
}

void obtenir_la_base_prototype(cycles *liste,int nb_cycles, graphemol m, struct graphe_dr* g_dr_i)
{
	
	/*
	 *  classer les cycles par taille de poids ( base minimale selon horton)
	 * */
	int i,j,l,k;

	int position[nb_cycles];
	int newposition[nb_cycles];
	for ( i = 0; i < nb_cycles;i++)
	{	
		//printf("\n cycle %d: ",i);
		/*for(j = 0; j < liste[i].nb_atomes;j++){
			printf("%d ", liste[i].sommets[j]);
		}*/
		position[i] = liste[i].nb_atomes;
		newposition[i]= 0;
	}
	//printf("\n");
	int max =0;
	for ( i = 0; i < nb_cycles;i++)
	{
		if(liste[i].nb_atomes > max)
			max = liste[i].nb_atomes;
	}
	int min = 0;
	for(i = 1; i < nb_cycles;i++)
	{
		if( position[i] < position[min])
		{
			min = i;
		}
	}
	newposition[0] = min;
	position[min] = max * max;
	for(i = 1; i < nb_cycles;i++)
	{
		//trouver le premier non nul 
			for(j = 0; j < nb_cycles;j++)
			{
				if(position[j] != max *max )
				{
					min = j;
					break;
				}
			}
			for(j = min; j < nb_cycles;j++)
			{
				if( position[j] < position[min])
				{
					min = j;
				}

			}
			newposition[i ] = min;
			position[min] = max * max;

	}
	
	/*
	 *  préparation des matrices pour élimination de Gauss
	 * */
	int **matrice = malloc(nb_cycles * sizeof(int *));
	if( matrice == NULL)
		probleme_memoire();

	for ( i = 0; i < nb_cycles;i++)
	{
		matrice[i] = malloc(m.nb_liaisons  * sizeof(int));
		if( matrice[i] == NULL)
			probleme_memoire();
	}
	int **matrice2 = malloc(nb_cycles * sizeof(int *));
	if( matrice2 == NULL)
		probleme_memoire();

	for ( i = 0; i < nb_cycles;i++)
	{
		matrice2[i] = malloc(m.nb_liaisons  * sizeof(int));
		if( matrice2[i] == NULL)
			probleme_memoire();
	}

	//initialisation de la matrice 
	for ( i = 0; i < nb_cycles;i++)
	{
		for ( j = 0; j < m.nb_liaisons;j++){
			matrice[i][j] = 0;
		}
	}

	int *chemin;
	int taille,pos ;
	int tab_base[nb_cycles];
	// tab_base indique quels cycles appartiennent à l'union des bases de cycles (pour l'instant aucun)
	for(i=0;i<nb_cycles;i++)
	{
		tab_base[i] = 0;
	}
	for ( i = 0; i < nb_cycles;i++)
	{
		// je travaille sur le cycle liste[newposition[i]]
		//printf("working on %d \n", liste[newposition[i]].id_cycle);
		taille = liste[newposition[i]].chemin1[0] +  liste[newposition[i]].chemin2[0] + 1;
		chemin = concatener_deux_chemins(liste[newposition[i]]);
		for( j = 0; j < taille ;j++ )
		{	
			pos = position_de_arete(chemin[j],chemin[j + 1],m);
			matrice[i][pos] = 1 ;
		}
		free(chemin);
	}

	//sauvegarde 
	for (i = 0; i < nb_cycles; i++)
	{
		for (j= 0; j < m.nb_liaisons; j++)
		{
			matrice2[i][j] = matrice[i][j];
		}
	}

	/*
	 * calcul des cycles relevants parmi l'ensemble de prototypes (élimination de Gauss)
	 * */
	
	// on réinitialise matrice
	int first;
	for(i=1;i < nb_cycles;i++)
	{
		for(j=0;j<m.nb_liaisons;j++)
		{
			matrice[i][j] = 0;
		}
	}

	int compte, compte_2, compte_k, compte_j;
	for(i=0;i < nb_cycles;i++)
	{
		// on teste le ième cycle
		compte = 0;
		for(j=0;j<m.nb_liaisons;j++)
		{
			//printf("matrice i j : %d\n", matrice2[i][j]);
			matrice[i][j] = matrice2[i][j];
			if(matrice[i][j] == 1) compte ++;
		}
		// Gauss sur B< et C
		k = 0;
		while(k < nb_cycles)
		{
			first = -1;
			compte_k = 0;
			for(j=0;j<m.nb_liaisons;j++)
			{
					if(matrice[k][j] == 1)
					{
						if(first == -1) first = j;
					}
					if(matrice2[k][j] == 1)
					{
						compte_k ++;
					}
			}
			//printf("compte : %d, compte_k :%d\n", compte, compte_k);
			if(compte_k >= compte && k != i)
			{
				k++;
				continue;
			}
			//printf("k : %d\n", k);
			
			if(first != -1) 
			{
				for ( j = k+1; j < nb_cycles; j++)
				{
					compte_j = compte_1_ligne(matrice2, m.nb_liaisons, j);
					if( matrice[j][first] == 1 && (compte_j < compte || j == i))
					{
						//printf("%d %d\n", j, k);
						
						for ( l = 0; l < m.nb_liaisons;l++)
						{  
							matrice[j][l] = fonction_xor(matrice[j][l],matrice[k][l]);
						}
					}
				}
			}
			k++;
		}
		compte_2 = 0;
		for(j=0;j<m.nb_liaisons;j++)
		{
			//printf("%d ", matrice[i][j]);
			if(matrice[i][j] == 1) compte_2 ++;
		}
				
		if(compte_2 != 0) // on garde le cycle C
		{
				//printf("gardes\n");
				tab_base[i] = 1;
				for(j=0;j<m.nb_liaisons;j++)
				{
					matrice[i][j] = matrice2[i][j];
				}
		}

		for(k=0;k<i;k++)
		{
			if(tab_base[k] == 1)
			{
				for(j=0;j<m.nb_liaisons;j++)
				{
					matrice[k][j] = matrice2[k][j];
				}
			}
		}
	}
		

	// on regarde quels cycles garder
	int somme;
	int base = 0;
	for ( i = 0; i < nb_cycles; i++)
		tab_base[i] = 0;
	for (i = 0; i < nb_cycles; i++)
	{
		somme = 0;
		for (j= 0; j < m.nb_liaisons; j++)
		{
			somme += matrice[i][j];
		}
		if( somme != 0)
		{
			base++;
			tab_base[i] = 1;
			//printf("gardes: %d\n", i+1);
		}
	}

	/* 
	 *  crée les cycles de chaque famille des prototypes relevants
	 * */
	int* chemin1, *chemin2;
	int s, nb_chemins_courants_1, nb_chemins_courants_2;
	int** chemins_p;
	int** chemins_q;

	cycles c;
	int nb_new_cycles = 0;
	
	for(i=0;i<nb_cycles;i++)
	{
		if(tab_base[i] == 1)
		{
			nb_chemins_courants_1 = 0;
			nb_chemins_courants_2 = 0;
			chemin1 = NULL;
			chemin2 = NULL;
			
			chemins_p = NULL;
			chemins_q = NULL;

			s =  liste[newposition[i]].sommet;
			//printf("s : %d, a1 : %d, a2 : %d\n", s, liste[newposition[i]].arete.A1,liste[newposition[i]].arete.A2);
		
			chemins_p = list_path(liste[newposition[i]].arete.A1, s, chemin1, g_dr_i[position_graphemol(m, s)], &nb_chemins_courants_1, liste[newposition[i]].chemin1[0]+1, chemins_p, liste[newposition[i]].arete.A1);
			chemins_q = list_path(liste[newposition[i]].arete.A2, s, chemin2, g_dr_i[position_graphemol(m, s)], &nb_chemins_courants_2, liste[newposition[i]].chemin2[0]+1, chemins_q, liste[newposition[i]].arete.A2);

			// on calcule les nouveaux cycles potentiels à partir des ensembles de chemins dans chemins_p et chemins_q
			for(j=0;j<nb_chemins_courants_1;j++)
			{
				
				chemin1 = malloc(chemins_p[j][0]*sizeof(int));
				chemin1[0] = chemins_p[j][0];
				for(l=2;l<chemins_p[j][0]+1;l++)
				{
					chemin1[l-1] = chemins_p[j][l];
				}
				
				for(k=0;k<nb_chemins_courants_2;k++)
				{
					chemin2 = malloc(chemins_q[k][0]*sizeof(int));
					chemin2[0] = chemins_q[k][0];
					for(l=2;l<chemins_q[k][0]+1;l++)
					{
						chemin2[l-1] = chemins_q[k][l];
					}
					
					if(chemin_independant(chemin1,chemin2,liste[newposition[i]].arete))
					{
						c = creer_un_cycle(m,liste[newposition[i]].sommet,liste[newposition[i]].arete,chemin1,chemin2);
						c.pere = m.nb_connexe;
						if(!verification_ajout_cycle(liste,nb_cycles+nb_new_cycles,c,m))
						{
							liste = ajouter_un_cycle(liste,nb_cycles+nb_new_cycles,c,m);
							nb_new_cycles++;
						}
						liberer_un_cycle(c);
					}
					free(chemin2);
				}
				free(chemin1);
			}

			for(j = 0; j< nb_chemins_courants_1;j++)
			{
				free(chemins_p[j]);
			}
			for(j = 0; j< nb_chemins_courants_2;j++)
			{
				free(chemins_q[j]);
			}
			free(chemins_p);
			free(chemins_q);
		}
	}
	
	
	// ajout des premiers cycles trouvés (les prototypes)
	for(i = 0; i < nb_cycles;i++)
	{
		if(tab_base[i] == 1)
		{
			//printf("ajouter cycle base: %d\n", i+1);
			ajouter_cycle_base(liste[newposition[i]],m);
		}
	}
	// ajout des cycles obtenus par famille de prototypes
	for(i=nb_cycles;i< nb_cycles+nb_new_cycles;i++)
	{
		//affichage_cycle(liste[i]);
		ajouter_cycle_base(liste[i],m);
	}
	
	// désallocation
	for(i = 0; i < nb_cycles+nb_new_cycles;i++)
	{
		liberer_un_cycle(liste[i]);
	}
	free(liste);
	
	for(i = 0; i < nb_cycles;i++)
	{
		free(matrice[i]);
		free(matrice2[i]);
	}
		
	free(matrice);
	free(matrice2);
	
}


int position_graphemol_arete(graphemol g, int sommet1,int sommet2)
{

	int pos = -1,i;
	for( i = 0; i < g.nb_liaisons;i++)
	{
		if((g.liste_liaisons[i].A1 == sommet1 && g.liste_liaisons[i].A2 == sommet2)||(g.liste_liaisons[i].A1 == sommet2 && g.liste_liaisons[i].A2 == sommet1))
		{
			pos = i;
			break;
		}
	}
	if(pos == -1)
	{
		fprintf(stdout, "liaison non trouve dans la molecule \n" );
		exit(555);
	}
	return pos;
}
void arete_dans_cycle(struct molecule m)
{
	graphemol g = conversion_mol_graphe(m);
	arete_cycle = malloc(g.nb_liaisons * sizeof(int));
	int i,j;
	for ( i = 0; i < g.nb_liaisons;i++ )
		arete_cycle[i] = 0;
	int *chemin,taille;
	for(i = 0; i < taille_base;i++)
	{
		taille = labase[i].chemin1[0] +  labase[i].chemin2[0] + 1;
		
		chemin = concatener_deux_chemins(labase[i]);
		for( j = 0; j < taille ; j++)
		{
			arete_cycle[position_graphemol_arete(g,chemin[j],chemin[j+1])]++;
		}

		free(chemin);
	}

	liberer_graphemol(g);
}

void arete_dans_cycle_liste(struct molecule m)
{
	graphemol g = conversion_mol_graphe(m);
	arete_liste = malloc(g.nb_liaisons *sizeof(int *));
	if(arete_liste == 0)
		probleme_memoire();
	int i;
	for ( i = 0; i < g.nb_liaisons;i++)
	{
		arete_liste[i] = malloc(arete_cycle[i]* sizeof(int));
		if( arete_liste[i] ==  NULL)
			probleme_memoire();
	}
	int *chemin,taille,k,j,pos;
	for ( k = 0; k < g.nb_liaisons;k++)
	{
		pos =0;
		//printf("-- %d :", k);
		for(i = 0; i < taille_base;i++)
		{
			taille = labase[i].chemin1[0] +  labase[i].chemin2[0] + 1;
			
			chemin = concatener_deux_chemins(labase[i]);
			for( j = 0; j < taille ; j++)
			{
				if((g.liste_liaisons[k].A1 == chemin[j] && g.liste_liaisons[k].A2 == chemin[j+1])||(g.liste_liaisons[k].A1 == chemin[j+1] && g.liste_liaisons[k].A2 == chemin[j]))
				{
					arete_liste[position_graphemol_arete(g,chemin[j],chemin[j+1])][pos] = i;
					pos++;
				}
			}

			free(chemin);
		}

	}

	//affichage
	liberer_graphemol(g);

}	
int sommets_commun ( cycles a, cycles b)
{
	int i,j,somme = 0;

	for (i = 0; i < a.nb_atomes; i++)
	{
		for (j = 0; j < b.nb_atomes; j++)
		{
			if(a.sommets[i] == b.sommets[j])
				somme++;
		}
	}


	return somme;
}
void nouvelle_arete(ARETE a)
{
	if(nb_arete_base == 0)
	{
		base_aretes= malloc((nb_arete_base+1)* sizeof(ARETE));
	}
	else
	{
		base_aretes= realloc(base_aretes,(nb_arete_base+1)* sizeof(ARETE));
	}
	if( base_aretes == NULL)
		probleme_memoire();
	base_aretes[nb_arete_base].id1=a.id1;
	base_aretes[nb_arete_base].id2=a.id2;
	base_aretes[nb_arete_base].type=a.type;
	base_aretes[nb_arete_base].poids=a.poids;

	nb_arete_base++;
}

int distance_inter_magma(cycles a , cycles b , struct molecule m)
{
	graphemol h = conversion_mol_graphe(m);
	int dist = h.nb_atomes * h.nb_atomes;
	int somme;
	struct liste_voisins * v = construction_voisinage_graphemol(h);
	if(!existe_chaine_graphemol(position_graphemol(h,a.sommet),position_graphemol(h,b.sommet),h,v))
		dist = -1;
	else
	{	
		int i,j,l,*pcc; 
		for( i = 0; i < a.nb_atomes;i++)
		{
			for( j = 0; j < b.nb_atomes;j++)
			{
				pcc = plus_court_chemin(a.sommets[i],b.sommets[j],h);
				if(pcc[0] < dist)
				{
					somme= 0;
					if(pcc[0] > 2)
					{
						somme += arete_cycle[position_graphemol_arete(h,a.sommets[i],pcc[pcc[0]-1])];
						somme += arete_cycle[position_graphemol_arete(h,pcc[1],b.sommets[j])];
						//printf("%d %d \n",a.sommets[i],pcc[pcc[0]-1]);
						//printf("%d %d --- > %d\n",a.sommets[i],pcc[1],arete_cycle[position_graphemol_arete(h,a.sommets[i],pcc[pcc[0]-1])] );
					}
					for( l = 1; l < pcc[0] - 1 ;l++)
						somme += arete_cycle[position_graphemol_arete(h,pcc[l],pcc[l+1])];
					if(somme == 0)
						dist = pcc[0];
				}
					
				free(pcc);
			}
		}
	}
	liberer_memoire_voisins(v,h);
	liberer_graphemol(h);
	return dist;
}


int nombre_de_cycles(int sommet1,int *chemin,int sommet2,graphemol g)
{
	int taille = 0;
	int tableau[taille_base];
	int i,r;
	for( i = 0; i < taille_base;i++)
		tableau[i] = 0;
	if(chemin[0] == 1)
	{
		taille = arete_cycle[position_de_arete(sommet1,sommet2,g)];
	}
	else
	{
		if(arete_cycle[position_de_arete(sommet1,chemin[chemin[0]-1],g)] > 0)
		{
			for( i = 0; i< arete_cycle[position_de_arete(sommet1,chemin[chemin[0]-1],g)] ;i++)
			{
				tableau[arete_liste[position_de_arete(sommet1,chemin[chemin[0]-1],g)][i]] = 1;
			}
		}
		//milieu a faire 
		if(chemin[0]>=3)
		{
			for( r = 1; r <chemin[0]-1;r++)
			{
				if(arete_cycle[position_de_arete(chemin[r],chemin[r+1],g)] > 0)
				{
					for( i = 0; i< arete_cycle[position_de_arete(chemin[r],chemin[r+1],g)] ;i++)
					{
						tableau[arete_liste[position_de_arete(chemin[r],chemin[r+1],g)][i]] = 1;
					}
				}

			}

		}
		if(chemin[0] >= 2)
		{
			if(arete_cycle[position_de_arete(sommet2,chemin[1],g)] > 0)
			{
				for( i = 0; i< arete_cycle[position_de_arete(sommet2,chemin[1],g)] ;i++)
				{
					tableau[arete_liste[position_de_arete(sommet2,chemin[1],g)][i]] = 1;
				}
			}
		}
		for (i = 0; i < taille_base;i++)
			taille += tableau[i];
	}
	return taille;
}


int verification_LC(cycles a, cycles b,struct molecule m,int dist)
{
	
	graphemol g = conversion_mol_graphe(m);
	int res = taille_base * taille_base;
	int i,j;
	int *chemin;
	int valeur,distance = g.nb_liaisons;
	for ( i = 0; i < a.nb_atomes;i++)
	{
		for(j = 0; j < b.nb_atomes;j++)
		{
			chemin = plus_court_chemin(a.sommets[i],b.sommets[j],g);
			
			if(a.sommets[i] != b.sommets[j])
			{
				valeur = nombre_de_cycles(a.sommets[i],chemin,b.sommets[j],g);
				if(valeur == 0 && distance >= chemin[0])
				{
					res = valeur;
					distance  = chemin[0];
				}
			}
			free(chemin);
		}
	}
	liberer_graphemol(g);
	if(res == taille_base * taille_base || distance != dist)
		res  = -1;
	return res;
}

// renvoie 1 si les chemins1 et 2 n'ont aucun sommet en commun et 0 sinon
int intersection_vide_chemins(int* chemin1, int* chemin2, int ext1, int ext2)
{
	int i,j;
	for(i=1;i<chemin1[0];i++)
	{
		for(j=1;j<chemin2[0];j++)
		{
			if(chemin1[i] == chemin2[j]) return 0;
		}
	}
	
	for(i=1;i<chemin1[0];i++)
	{
		if(ext2 == chemin1[i]) return 0;
	}
	
	for(i=1;i<chemin2[0];i++)
	{
		if(ext1 == chemin2[i]) return 0;
	}
	
	return 1;
}

void libere_graphe_dr(struct graphe_dr* g_dr, int nb_atomes)
{
	int i;
	for(i=0;i<nb_atomes;i++)
	{
		free(g_dr[i].liste_sommets);
		free(g_dr[i].liste_arcs);
	}
	free(g_dr);
	
}

int est_dans_liste_sommets(SOMMET_VR* liste_s, int nb_sommet_vr, int x)
{
	int i;
	for(i=0;i<nb_sommet_vr;i++)
	{
		if(liste_s[i].id == x) return 1;
	}
	return 0;
}

void trouver_prototypes_cycles_vismara(graphemol t)
{	

	int* chemin;
	int i,j,k,l,m,pos_k, pos_l, taille, nb_sommets_vr, nb_sommets_vr_temp;
	struct liaison li;
	int* liste_sommets; // ensemble de sommets Vi
	int** liste_chemins;
	int* S; // ensemble de sommets
	
	struct liste_voisins *v = construction_voisinage_graphemol(t);
	
	cycles *liste = NULL;
	cycles c;
	int nb_cycles  = 0; 
	
	struct graphe_dr* g_dr = malloc(t.nb_atomes*sizeof(struct graphe_dr));
	ARC* liste_arcs = NULL;
	SOMMET_VR* liste_sommets_vr = NULL;
	
	int nb_voisins_i = 0;
	for(i=0;i<t.nb_atomes;i++)
	{
		//printf("atomes %d\n", t.liste_atomes[i]);
		liste_sommets = malloc(t.nb_atomes*sizeof(int));
		for(j=0;j<t.nb_atomes;j++)
		{
			liste_sommets[j] = 0;
		}
		liste_chemins = malloc(t.nb_atomes*sizeof(int*));
	
		// calcul des plus courts chemins de i à tous les sommets j si ces plus courts chemins ne contiennent que des sommets précédant i (dans l'ordre croissant)
		// et calcul de l'ensemble Vr des sommets j ayant cette propriété
		for(j=0;j<t.nb_atomes;j++)
		{
				liste_chemins[j] = plus_court_chemin_precedant(t.liste_atomes[i], t.liste_atomes[j], t);
				// stockage d'un plus court chemin dans liste_chemins[j] = [nb aretes du plus court chemin entre i et j, voisin de j dans le chemin, voisin du voisin de j, ..., voisin de i dans le chemin]

				if(liste_chemins[j] != NULL && liste_chemins[j][0] > 0)
				{
					liste_sommets[j] = 1;
				}
		}
		
		taille = 0;
		nb_sommets_vr = 0;
		for(j=0;j<t.nb_atomes;j++)
		{
			//printf("liste_atomes[j] : %d\n", t.liste_atomes[j]);
			if(liste_sommets[j])
			{
				nb_sommets_vr_temp = nb_sommets_vr;
				if(!est_dans_liste_sommets(liste_sommets_vr, nb_sommets_vr, t.liste_atomes[j]))
				{
					nb_sommets_vr++;
					for(k=0;k<v[i].nb_voisins;k++)
					{
						if(v[i].id_voisins[k] == t.liste_atomes[j]) nb_voisins_i ++;
					}
				}
				
				S = malloc(t.nb_atomes*sizeof(int)); // va contenir: ensemble des sommets appartenant à un plus court chemin de i à j
				for(k=0;k<t.nb_atomes;k++)
				{
					S[k] = -1;
				}
				m = 0;
				
				for(k=0;k<v[j].nb_voisins;k++)
				{
					pos_k = position_graphemol(t, v[j].id_voisins[k]);
					if(liste_sommets[pos_k]) // le voisin k du sommet j appartient à Vi
					{
						if(liste_chemins[pos_k][0] + 1 == liste_chemins[j][0])
						{ // le sommet k appartient à un plus court chemin de i à j
							S[m] = v[j].id_voisins[k];
							//printf("S[m] : %d\n", S[m]);
							if(!est_dans_liste_sommets(liste_sommets_vr, nb_sommets_vr_temp, S[m]))
							{
								nb_sommets_vr++;
								for(l=0;l<v[i].nb_voisins;l++)
								{
									if(v[i].id_voisins[l] == S[m]) nb_voisins_i ++;
								}
							}
							m++;
						}
						else if(liste_chemins[pos_k][0] + 1 != liste_chemins[j][0] && pos_k < j && intersection_vide_chemins(liste_chemins[j], liste_chemins[pos_k], t.liste_atomes[j], v[j].id_voisins[k]) && liste_chemins[pos_k][0] == liste_chemins[j][0])
						{ // on ajoute un cycle de taille impaire

							for(l=0;l<t.nb_liaisons;l++)
							{
								if(t.liste_liaisons[l].A1 == t.liste_atomes[j] && t.liste_liaisons[l].A2 == t.liste_atomes[pos_k]) 
								{
									li = t.liste_liaisons[l];
									c = creer_un_cycle(t,t.liste_atomes[i],li,liste_chemins[j],liste_chemins[pos_k]);
								}
								else if (t.liste_liaisons[l].A2 == t.liste_atomes[j] && t.liste_liaisons[l].A1 == t.liste_atomes[pos_k]) 
								{
									li = t.liste_liaisons[l];
									c = creer_un_cycle(t,t.liste_atomes[i],li,liste_chemins[pos_k],liste_chemins[j]);
								}
							}
							
							c.pere = t.nb_connexe;
							if(!verification_ajout_cycle(liste,nb_cycles, c,t))
							{
								liste = ajouter_un_cycle(liste,nb_cycles,c,t);
								nb_cycles++;
							}
							liberer_un_cycle(c);
						}
						
					}
				}
				
				// liste_sommets_vr et liste_arcs sont pour construction du graphe Dr (voir article Vismara) utile pour rechercher tous les cycles d'une famille 
				// graphe Dr stocke tous les plus courts chemins depuis un sommet r (ne contenant que des sommets de numéro plus petits que r)
				if(liste_arcs == NULL) liste_arcs = malloc((taille+m)*sizeof(ARC));
				else liste_arcs = realloc(liste_arcs, (taille+m)*sizeof(ARC));
				//printf("nb %d\n", nb_sommets_vr_temp);
				
				if(liste_sommets_vr == NULL) liste_sommets_vr = malloc(nb_sommets_vr*sizeof(SOMMET_VR));
				else liste_sommets_vr = realloc(liste_sommets_vr, nb_sommets_vr*sizeof(SOMMET_VR));
				
				//printf("nb sommets vr temp : %d\n", nb_sommets_vr_temp);
				if(!est_dans_liste_sommets(liste_sommets_vr, nb_sommets_vr_temp, t.liste_atomes[j])) 
				{
					liste_sommets_vr[nb_sommets_vr_temp].id = t.liste_atomes[j];
					nb_sommets_vr_temp ++;
				}
				
				for(l=0;l<m;l++)
				{
					if(!est_dans_liste_sommets(liste_sommets_vr, nb_sommets_vr_temp, S[l]))
					{ 
						liste_sommets_vr[nb_sommets_vr_temp].id = S[l];
						nb_sommets_vr_temp ++;
					}
					liste_arcs[taille+l].id1 = t.liste_atomes[j];
					liste_arcs[taille+l].id2 = S[l];
				}
				taille = taille + m;
				//printf("nb sommets : %d %d\n", nb_sommets_vr, nb_sommets_vr_temp);
				// entre toute paire de sommets de S, on ajoute un cycle de taille paire si les chemins depuis i sont disjoints
				k=0;
				while(S[k] != -1)
				{
					//printf("S[k] %d\n", S[k]);
					pos_k = position_graphemol(t, S[k]);
					l = k+1;
					while(S[l] != -1)
					{
						//printf("S[l] %d\n", S[l]);
						pos_l = position_graphemol(t, S[l]);
						if(intersection_vide_chemins(liste_chemins[pos_k], liste_chemins[pos_l], S[k], S[l]))
						{
							//printf("taille chemin : %d\n", liste_chemins[pos_k][0]);
							//printf("%d \n", S[k]);
							chemin = malloc((liste_chemins[pos_k][0]+1)*sizeof(int));
							
							chemin[0] = liste_chemins[pos_k][0]+1;
							chemin[1] = S[k];
							for(m=1;m<liste_chemins[pos_k][0];m++)
							{
								chemin[m+1] = liste_chemins[pos_k][m];
							}
							
							for(m=0;m<t.nb_liaisons;m++)
							{
								if(t.liste_liaisons[m].A1 == t.liste_atomes[j] && t.liste_liaisons[m].A2 == t.liste_atomes[pos_l]) 
								{
									li = t.liste_liaisons[m];
									c = creer_un_cycle(t,t.liste_atomes[i],li,chemin,liste_chemins[pos_l]);
								}
								else if (t.liste_liaisons[m].A2 == t.liste_atomes[j] && t.liste_liaisons[m].A1 == t.liste_atomes[pos_l])
								{
									li = t.liste_liaisons[m];
									c = creer_un_cycle(t,t.liste_atomes[i],li,liste_chemins[pos_l], chemin);
								}
							}
							
							c.pere = t.nb_connexe;
							if(!verification_ajout_cycle(liste,nb_cycles, c,t))
							{
								liste = ajouter_un_cycle(liste,nb_cycles,c,t);
								nb_cycles++;
							}
							liberer_un_cycle(c);
							free(chemin);
							
						} 
						l++;
					}
					k++;
				}
				free(S);
			}
		}
		// on rajoute le graphe Dr du sommet i courantà l'ensemble des graphes Dr (dans g_dr)
		nb_sommets_vr ++;
		if(liste_sommets_vr == NULL) liste_sommets_vr = malloc(nb_sommets_vr*sizeof(SOMMET_VR));
		else liste_sommets_vr = realloc(liste_sommets_vr, nb_sommets_vr*sizeof(SOMMET_VR));
		liste_sommets_vr[nb_sommets_vr-1].id = t.liste_atomes[i];
		
		m = taille;
		if(liste_arcs == NULL) liste_arcs = malloc(nb_voisins_i*sizeof(ARC));
		else liste_arcs = realloc(liste_arcs, nb_voisins_i*sizeof(ARC));

		//printf("nb sommets vr : %d\n", nb_sommets_vr);
		//printf("nb voisins i : %d\n", v[i].nb_voisins);
		for(j = 0;j<nb_sommets_vr;j++)
		{
			//printf("%d\n", liste_sommets_vr[j].id);
			for(k=0;k<v[i].nb_voisins;k++)
			{
				//printf("%d %d\n", j, k);
				if(v[i].id_voisins[k] == liste_sommets_vr[j].id)
				{
					liste_arcs[m].id1 = liste_sommets_vr[j].id;
					liste_arcs[m].id2 = t.liste_atomes[i];
					m ++;
				}
			}
		}
		//printf("nb sommets : %d, nb_arcs : %d\n", nb_sommets_vr, m);
		//printf("liste sommets : %p %d\n", liste_sommets_vr, liste_sommets_vr[0].id);
		g_dr[i].nb_sommets = nb_sommets_vr;
		g_dr[i].nb_arcs = m;
		g_dr[i].liste_sommets = (SOMMET_VR*) malloc(nb_sommets_vr*sizeof(SOMMET_VR));
		for(j = 0; j < nb_sommets_vr;j++)
		{
			g_dr[i].liste_sommets[j] = liste_sommets_vr[j];
		}
		g_dr[i].liste_arcs = malloc(m*sizeof(ARC));
		for(j = 0; j < m;j++)
		{
			g_dr[i].liste_arcs[j].id1 = liste_arcs[j].id1;
			g_dr[i].liste_arcs[j].id2 = liste_arcs[j].id2;
		}
		
		free(liste_sommets_vr);
		free(liste_arcs);
		liste_sommets_vr = NULL;
		liste_arcs = NULL;
		//printf("nb sommets : %d, nb_arcs : %d\n", g_dr[i].nb_sommets, m);
		//printf("liste sommets : %p %d\n", g_dr[i].liste_sommets, g_dr[i].liste_sommets[0].id);
		
		
		
		free(liste_sommets);
		for(j=0;j<t.nb_atomes;j++)
		{
			if(liste_chemins[j] != NULL) free(liste_chemins[j]);
		}
		free(liste_chemins);
	}
	
	// calcule l'union des bases de cycles
	obtenir_la_base_prototype(liste, nb_cycles, t, g_dr); 
	
	// désallocation
	liberer_memoire_voisins(v,t);
	libere_graphe_dr(g_dr, t.nb_atomes);

}

/* Prétraitement (élimination des sommets de degré 0 ou 1 et des isthmes et décomposition en composantes connexes)
 * Calcul de l'union des bases de cycles (Vismara 97)
*/
void elimination_feuilles(struct molecule m)
{
	
	graphemol g = conversion_mol_graphe(m);
	
	// enlever les sommets qui ne sont dans aucun cycle
	int *deja_elimine = malloc(g.nb_atomes * sizeof(int));
	int i;
	for (i = 0; i < g.nb_atomes; i++)
		deja_elimine[i] = 0;
	int *degre = calcul_degre_mol(g);
	g = modification_structure_mol(g,degre,deja_elimine);

	int sommets  = g.nb_atomes;
	for ( i = 0;  i < g.nb_atomes;i++)
	{
		//printf("-%d %d\n", i+1,deja_elimine[i]);
		if( deja_elimine[i] == 1)
			sommets--;
	}
	int l;
	
	if ( sommets > 0)
	{	
		/* 
		 * Recherche l'ensemble des composantes connexes du graphe
		 * (sans les atomes qui sont pas dans des cycles)
		 * */
		graphemol *liste_connexe = ensemble_connexe_graphemol(g, deja_elimine);

		//printf("connexe = %d\n",liste_connexe[0].nb_connexe);
		
		/* 
		 * Recherche tous les isthmes et les enlève
		 * */
		int nb_connexe = liste_connexe[0].nb_connexe;
		int nb_isthmes = 0;
		int *elimine;
		for (i = 0; i < nb_connexe; i++)
		{
			//printf("i = %d\n", i + 1);
			//affiche_graphemol(liste_connexe[i]);
			elimine = malloc(liste_connexe[i].nb_atomes * sizeof(int));
			for ( l = 0;  l < liste_connexe[i].nb_atomes;l++)
			{	
				elimine[l] = 0;
			}

			nb_isthmes += nombre_isthmes_graphemol(liste_connexe[i],elimine);
			free(elimine);
		}
		
		isthmes *liste_isthmes = malloc(nb_isthmes *sizeof(isthmes));
		if( liste_isthmes == NULL)
			probleme_memoire();

		liste_isthmes = retrouver_tous_isthmes(liste_connexe,liste_isthmes,deja_elimine);

		//for(i = 0; i< nb_isthmes;i++)
		//	printf("%d - %d %d\n",i+1,liste_isthmes[i].l.A1,liste_isthmes[i].l.A2 );
		graphemol h = enlever_tous_isthmes(g,deja_elimine,liste_isthmes,nb_isthmes);
		
		deja_elimine = realloc(deja_elimine,h.nb_atomes *sizeof(int));
		for (i = 0; i < h.nb_atomes; i++)
			deja_elimine[i] = 0;
		
		/* 
		 * Recalcule les composantes connexes après avoir enlevé les isthmes
		 * */
		 
		graphemol *t = ensemble_connexe_graphemol(h,deja_elimine);
		//printf("nb total composants bas niveau %d at = %d\n",t[0].nb_connexe,t[0].nb_atomes);
		//printf("connexe de t = %d\n",t[0].nb_connexe);
		//affiche_graphemol(t[0]);
		
		/*
		 * Trouve l'union des bases de cycles (Vismara, 1997)
		 * */
		int p,j;
		for (i = 0; i < t[0].nb_connexe; i++)
		{
			p = t[i].liste_atomes[0];
			for( j = 0; j < nb_connexe;j++)
			{
				if(sommet_dans_basniveau(liste_connexe[j],p))
				{
					//printf("p: %d est dans %d\n",p,j+1 );
					t[i].pere = j;
					break;
				}
			}
			trouver_prototypes_cycles_vismara(t[i]);
		}
		
		//trouver le composant bas niveauu de chaque cycle de la base
		//printf("taille de la base : %d \n",taille_base);
		for( i = 0; i < taille_base; i++)
		{
			p = labase[i].sommets[0];
			for( j = 0; j < t[0].nb_connexe;j++)
			{
				if(sommet_dans_basniveau(t[j],p))
				{
					//printf("p: %d est dans %d et papy %d\n",p,j+1,t[j].pere );
					labase[i].pere = j;
					break;
				}
			}
		}
		
		//creation des aretes
		arete_dans_cycle(m);
		arete_dans_cycle_liste(m);
		int dist,dist2;
		ARETE a;
		for( i = 0; i < taille_base - 1; i++)
		{
			for( j = i+1; j < taille_base; j++)
			{
				if(labase[i].pere == labase[j].pere) //arete de type 1 ou 2 meme magma
				{
					if(sommets_commun(labase[i],labase[j]) >= 1) // au moins un sommet en commun ---- bleu
					{
						a.id1 = i;
						a.id2 = j;
						a.type = 1;
						a.poids = sommets_commun(labase[i],labase[j]) -1;
						nouvelle_arete(a);
					}
				}
				else //arete de type 3 magma different 
				{// existe til un chemin de c1 vers c2 
					dist = distance_inter_magma(labase[i],labase[j],m); //----- vert
					if( dist !=-1)
					{//il existe un chemin entre ces deux cycles dans la molecule 
						dist2 = verification_LC(labase[i],labase[j],m,dist);
						if( dist2 != -1)
						{
							a.id1 = i;
							a.id2 = j;
							a.type = 3;
							a.poids = dist;
							nouvelle_arete(a);
						}
					}
				}		
			}
		}


		//liberation de la memoire
		for (i = 0; i < nb_connexe; i++)
			liberer_graphemol(liste_connexe[i]);
		nb_connexe = t[0].nb_connexe;
		for ( i = 0; i < nb_connexe; i++)
			liberer_graphemol(t[i]);	
		free(liste_connexe);
		free(liste_isthmes);
		liberer_graphemol(h);
		free(t);
		//printf("TAILLE DE LA BASE : %D ",taille_base);
		
		graphemol sd = conversion_mol_graphe(m);
		if( sd.nb_liaisons > 0)
		{
			for( i = 0; i < sd.nb_liaisons;i++)
				free(arete_liste[i]);
			free(arete_liste);
			free(arete_cycle);
		}
		liberer_graphemol(sd);

	}

	
	
	
	free(deja_elimine);
		
	liberer_graphemol(g);
	

	
}


ARETE copier_arete(ARETE a , ARETE b)
{
	a.id1 = b.id1;
	a.id2 = b.id2;
	a.type = b.type;
	a.poids = b.poids;

	return a;
}

GRAPHE_CYCLE construction_graphe_cycles(struct molecule m)
{
	nb_arete_base=0;
	taille_base=0;
	
	elimination_feuilles(m); // la recherche de cycles se fait dans cette fonction
	GRAPHE_CYCLE c;
	
	c.nb_sommets = taille_base;
	c.nb_aretes = nb_arete_base;
	int i;
	c.liste_sommets = NULL;
	c.liste_aretes  = NULL;
	c.liste_sommets = malloc( c.nb_sommets * sizeof(SOMMET));
	c.liste_aretes  = malloc(c.nb_aretes * sizeof(ARETE));
	//printf("nb sommets %d et aretes %d\n", c.nb_sommets,c.nb_aretes);
	for( i = 0; i < c.nb_sommets;i++)
	{
		//printf("%d id en cours %d\n",labase[i].id_cycle,labase[i].nb_atomes );
		c.liste_sommets[i].id = labase[i].id_cycle;
		c.liste_sommets[i].poids = labase[i].nb_atomes;
		c.liste_sommets[i].type = 1;
		c.liste_sommets[i].poids_bouboule = NULL;

	}
	for( i = 0; i < c.nb_aretes;i++)
	{
		c.liste_aretes[i] = copier_arete(c.liste_aretes[i],base_aretes[i]);
	}

	if(nb_arete_base > 0)
	{
		free(base_aretes);
	}

	return c;
}

void liberer_graphe_cycles( GRAPHE_CYCLE c)
{
	//int i ; 
	if(c.liste_sommets !=NULL)
		free(c.liste_sommets);
	if(c.liste_aretes !=NULL)
	free(c.liste_aretes);
	
}

int min( int a , int b)
{
	if ( a < b )
		return a;
	return b;
}
