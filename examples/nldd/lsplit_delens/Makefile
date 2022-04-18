#### Macros ####
exename=exe

## Complier and options ##
FC = ifort
FLAG = -O3 -ip -fpp
DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) #$(DBAG)

## Directories ##
libroot = ../../../
MyLib = ${libroot}/
LAPACK = ${libroot}/pub/LAPACK95

MODD = -I$(LAPACK)/mod -I$(MyLib)/mod
LIBD = -L$(LAPACK)/lib -L$(MyLib)/lib
LINK = -lnldd -llinalg -lutils -llapack95 -llapack -lrefblas
OPTION = $(MODD) $(LIBD) $(LINK)

OBJ = main.o

all: $(exename)

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

$(exename): $(OBJ)
	$(FC) $(FLAGS) $(OBJ) $(OPTION) -o $@

clean:
	rm -f *.o* *.e* *.mod *.d *.pc *.obj core* *.il

depend:
	makedepend -- $(F90FLAGS) -- $(OBJFILES)

