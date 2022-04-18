#////////////////////////////////////////////////////#
# * Makefile
#////////////////////////////////////////////////////#

## libfile ##
target=liblinalg.a

# complier and options
FC = ifort
FLAG = -O3 -fpp -fPIC
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

# directories
libroot = ../../
MyLib = ${libroot}
LAPACK = ${libroot}/pub/LAPACK95/

# compile links
MODD = -I$(MyLib)/mod -I$(LAPACK)/mod
LIBD = -L$(LAPACK)/lib
LINK = -lutils -llapack95 -llapack -lrefblas
OPTION = $(MODD) $(LIBD) $(LINK)

# files to be compiled
OBJ = myblas.o mylapack.o mylapack95.o

$(target): $(OBJ)
	ar crs $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	mv *.mod ${MyLib}/mod; mv $(target) ${MyLib}/lib/

clean:
	rm -f *.o* *.mod *.a

# file dependencies

