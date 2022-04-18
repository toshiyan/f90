#////////////////////////////////////////////////////#
# * Makefile
#////////////////////////////////////////////////////#

## libfile ##
target = libanaflat.a

## Complier and options ##
FC = ifort
FLAG = -O3 -fpp -fPIC
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

## Directories ##
MyLib = ../../

#FFTW3 = /global/common/cori/software/fftw/3.3.4/hsw/intel/
#MODD = -I$(FFTW3)/include -I$(MyLib)/mod
#LIBD = -L$(cfitsio) -L$(FFTW3)/lib -L$(MyLib)/lib
cfitsio = ${MyLib}/pub/cfitsio
FFTW3   = ${MyLib}/pub/FFTW/
MODD = -I$(FFTW3)/api -I$(MyLib)/mod
LIBD = -L$(cfitsio) -L$(FFTW3) -L$(MyLib)/lib
LINK = -lutils -lcfitsio -lfftw3

OPTION = $(MODD) $(LIBD) $(LINK)

OBJ = anaflat.o myfftw.o mycfitsio.o recflat.o rotrecflat.o taurecflat.o delens.o inp.o remapping.o bispecflat.o

$(target): $(OBJ)
	ar cr $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	mv *.mod ../../mod; mv $(target) ../../lib/ 

clean:
	rm -f *.o* *.mod *.a

myfftw.o : anaflat.o
mycfitsio.o : anaflat.o
recflat.o : myfftw.o
inp.o : myfftw.o
remapping.o : myfftw.o
bispflat.o : anaflat.o myfftw.o

