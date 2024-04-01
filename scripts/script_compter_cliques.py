import csv

with open("results/results_cliques_type_compte_reduit.csv", "w") as fic:
	csvwriter = csv.writer(fic)
	csvwriter.writerow(["Name", "Classe", "Bouboule", "Coin bouboule", "Coin pas ouvert", "Coin ouvert mais DL aux jonctions", "Coin ouvert"])
	with open("results/results_cliques_type_reduit.csv", "r") as fichier :
		csvreader = csv.reader(fichier)
		
		for row in csvreader :
			i = 2
			c_bouboule = 0
			c_coin_bouboule = 0
			c_coin_pas_ouvert = 0
			c_coin_ouvert = 0
			c_coin_ouvert_mais_dl = 0
			while row[i] != " " and row[i] != "" :
				if int(row[i]) == 0 :
					c_bouboule += 1
				elif int(row[i]) == 1 :
					c_coin_bouboule += 1
				elif int(row[i]) == 2 :
					c_coin_pas_ouvert += 1
				elif int(row[i]) == 3 : 
					c_coin_ouvert += 1
				elif int(row[i]) == 4 :
					c_coin_ouvert_mais_dl += 1
				i += 1
			csvwriter.writerow([row[0], row[1], c_bouboule, c_coin_bouboule, c_coin_pas_ouvert, c_coin_ouvert_mais_dl, c_coin_ouvert]) 
