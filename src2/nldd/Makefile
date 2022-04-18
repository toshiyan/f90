#////////////////////////////////////////////////////#
# * Makefile
#////////////////////////////////////////////////////#

## libfile ##
target=libnldd.a

## Complier and options ##
FC = ifort
FLAG = -O3 -fpp -fPIC
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

## Directories ##
libroot = ../../
MyLib = ${libroot}
LAPACK = ${libroot}/pub/LAPACK95/

MODD = -I$(LAPACK)/mod -I$(MyLib)/mod
LIBD = -L$(LAPACK)/lib -L$(MyLib)/lib
LINK = -llinalg -lutils -llapack95 -llapack -lrefblas

OPTION = $(MODD) $(LIBD) $(LINK)

OBJ = bb_flat.o nldd_kernel.o nldd_lens.o nldd_delens.o nldd_rot.o nldd_tau.o nldd_src.o nldd_interface.o nldd_flat.o #nldd_extkernel.o #nldd_ylm.o

$(target): $(OBJ)
	ar cr $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	mv *.mod ${MyLib}/mod; mv $(target) ${MyLib}/lib/

clean:
	rm -f *.o*

nldd_interface.o: nldd_kernel.o
nldd_ylm.o: nldd_interface.o


