## Exe file ##
exename=exe

## Complier and options ##
FC = ifort
FLAG = -O3 -fpp
FLAGS = $(FLAG) $(DBAG)

## Directories ##
libroot = ../../../
MyLib = ${libroot}
FFTW3 = ${libroot}/pub/FFTW

MODD = -I$(MyLib)/mod -I$(FFTW3)/api
LIBD = -L$(MyLib)/lib -L$(FFTW3)
LINK = -lanaflat -lutils -lfftw3

OPTION = $(MODD) $(LIBD) $(LINK)

OBJ = main.o

$(exename): $(OBJ)
	$(FC) $(FLAGS) $(OBJ) $(OPTION) -o $@

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

clean:
	rm -f *.o* *.mod

