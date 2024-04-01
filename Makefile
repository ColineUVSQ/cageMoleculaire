CFLAGS=-g -Wall
CC = gcc -Wall -Wextra
CXX = g++ -Wall -Wextra


run_cage: analyse_cage
	#valgrind -v --leak-check=full --show-leak-kinds=all ./analyse_cage
	./analyse_cage
	#gdb ./analyse_cage

analyse_cage: analyse_cage.o utils_cage_moleculaire.o lecture_molecule_sdf.o graphe_cycles.o
	gcc ${CFLAGS} -o $@ $^

analyse_cage.o: analyse_cage.c analyse_cage.h
	gcc ${CFLAGS} -c analyse_cage.c

utils_cage_moleculaire.o: utils_cage_moleculaire.c utils_cage_moleculaire.h structure.h
	gcc ${CFLAGS} -c utils_cage_moleculaire.c

graphe_cycles.o: graphe_cycles.c graphe_cycles.h
	gcc ${CFLAGS} -c graphe_cycles.c
	
lecture_molecule_sdf.o: lecture_molecule_sdf.c lecture_molecule_sdf.h 
	gcc ${CFLAGS} -c lecture_molecule_sdf.c


clean: 
	rm -f analyse_cage
	rm -f sortie
	rm -f *o

clean_results:
	rm -f -r results/
	rm -f -r data/smi_files_reduit/
	rm -f -r data/dot_files_reduit/
	rm -f -r data/png_files_reduit/

