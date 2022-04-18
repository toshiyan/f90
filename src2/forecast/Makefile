#////////////////////////////////////////////////////#
# * Makefile
#////////////////////////////////////////////////////#

## Exe file ##
target=libforecast.a

## Complier and options ##
FC = ifort
FLAG = -O3 -fpp -fPIC
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

## Directories ##
libroot = ../../../
MyLib = ${libroot}/mylib
LAPACK = ${libroot}/LAPACK95

MODD = -I$(LAPACK)/mod -I$(MyLib)/mod
LIBD = -L$(LAPACK)/lib -L$(MyLib)/lib
LINK = -llinalg -lutils -llapack95 -llapack -lrefblas

OPTION = $(MODD) $(LIBD) $(LINK)

OBJ = cosmoparams.o forecast.o results.o snr.o pca.o

$(target): $(OBJ)
	ar cr $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	mv *.mod ${MyLib}/mod; mv $(target) ${MyLib}/lib/ 

clean:
	rm -f *.o*

