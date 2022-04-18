## Exe file ##
exename=exe

## Complier and options ##
FC = ifort
FLAG = -O3 -fpp
FLAGS = $(FLAG) $(DBAG)

## Directories ##
libroot = ../../../
MyLib = ${libroot}

MODD = -I$(MyLib)/mod
LIBD = -L$(MyLib)/lib
LINK = -lutils

OPTION = $(MODD) $(LIBD) $(LINK)

OBJ = main.o

$(exename): $(OBJ)
	$(FC) $(FLAGS) $(OBJ) $(OPTION) -o $@

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

clean:
	rm -f *.o* *.mod

