FC=gfortran
# Liste d'options pour la phase de debug: à commenter en
# phase d'exploitation
FFLAGS=-O0 -g -Wall -ffpe-trap=invalid,overflow -fcheck=all
# Liste d'options d'optimisation du code  pour la phase
# d'exploitation': à commenter en phase de debug
# FFLGAS=-O3 -march=native
EXEC=run

all : $(EXEC)

# Edition de lien: construit l'exécutable à partie de tous les .o
# On met les dépendances (les .o) dans l'ordre dans lequel on veut
# qu'elles soient construites.
$(EXEC) : mod_precision.o mod_maillage.o mod_sortie.o main.o 
	$(FC) -o $@ $^ $(FFLAGS)

# Compilation : construire les .o à partir des .h
# attention, on ne contrôle pas la présence des .mod, il faut
# donc s'assurer de construire les .o dans le bon ordre.
%.o : %.f90
	$(FC) -c $< $(FFLAGS)

clean :
	rm -f *.o *.vtk *.mod $(EXEC)