## Exe file ##
exename=exe

## Complier and options ##
FC = ifort
FLAG = -O3 -fpp
FLAGS = $(FLAG) $(DBAG)

## Directories ##
MyLib = ../../../
FFTW3 = ${MyLib}/pub/FFTW
#FFTW3 = /global/common/cori/software/fftw/3.3.4/hsw/intel/
#MODD = -I$(MyLib)/mod -I$(FFTW3)/include
#LIBD = -L$(MyLib)/lib -L$(FFTW3)/lib
MODD = -I$(MyLib)/mod -I$(FFTW3)/api
LIBD = -L$(MyLib)/lib -L$(FFTW3)
LINK = -lanaflat -lutils -lfftw3

OPTION = $(MODD) $(LIBD) $(LINK)

OBJ = test.o

$(exename): $(OBJ)
	$(FC) $(FLAGS) $(OBJ) $(OPTION) -o $@

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

clean:
	rm -f *.o* *.mod

