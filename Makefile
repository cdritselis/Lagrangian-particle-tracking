sph: sph.o
	gfortran -o sph sph.o
sph.o: sph.F90
	gfortran -c -std=f2008 -Wall -g sph.F90
