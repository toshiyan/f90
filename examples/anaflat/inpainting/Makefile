## Exe file ##
exename=exe

## Complier and options ##
FC = ifort
FLAG = -O3 -ip -fpp
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

## Directories ##
MyLib = ../../../
healpix = ${MyLib}/pub/Healpix
LAPACK = ${MyLib}/pub/LAPACK95/lib
FFTW3 = ${MyLib}/pub/FFTW
cfitsio = ${MyLib}/pub/cfitsio

MODD = -I$(LAPACK) -I$(MyLib)/mod -I$(FFTW3)/api
LIBD = -L$(LAPACK) -L$(MyLib)/lib -L$(cfitsio) -L$(FFTW3)
LINK = -lanaflat -lutils -llapack95 -llapack -lrefblas -lfftw3 -lcfitsio

OPTION = $(MODD) $(LIBD) $(LINK)

OBJ = main.o

$(exename): $(OBJ)
	$(FC) $(FLAGS) $(OBJ) $(OPTION) -o $@

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

clean:
	rm -f *.o *.mod

