## Executable file ##
exename=exe

FC = ifort
FLAG = -O3 -ip -fpp -error-limit 5 -heap-arrays
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

libroot = ../../../
MyLib = ${libroot}
LAPACK = ${libroot}/pub/LAPACK95/

MODD = -I$(LAPACK)/mod -I$(MyLib)/mod
LIBD = -L$(MyLib)/lib -L$(LAPACK)/lib
LINK = -lnldd -llinalg -lutils -llapack95 -llapack -lrefblas

OBJ = main.o

all: $(exename)

$(exename): $(OBJ)
	$(FC) $(FLAGS) $(OBJ) $(LIBD) $(LINK) -o $@

%.o: %.f90
	$(FC) $(FLAGS) -c $(MODD) $*.f90

clean: 
	rm -f *.o* *.mod

