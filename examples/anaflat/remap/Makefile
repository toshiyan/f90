#### Macros ####
## Exe file ##
exename=exe

## Complier and options ##
FC = ifort
FLAG = -O3 -ip -fpp
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

## Directories ##
libroot = ../../../
MyLib = ${libroot}
FFTW3 = ${libroot}/pub/FFTW/

MODD = -I$(MyLib)/mod -I$(FFTW3)/api
LIBD = -L$(MyLib)/lib -L$(FFTW3)

LINK = -lanaflat -lutils -lfftw3
OPTION = $(MODD) $(LIBD) $(LINK)

OBJ = test.o

all: $(exename)

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

$(exename): $(OBJ)
	$(FC) $(FLAGS) $(OBJ) $(OPTION) -o $@

clean:
	rm -f *.o* *.e* *.mod *.d *.pc *.obj core* *.il

depend:
	makedepend -- $(F90FLAGS) -- $(OBJFILES)

