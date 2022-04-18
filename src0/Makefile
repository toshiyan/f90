#////////////////////////////////////////////////////#
# * Makefile
#////////////////////////////////////////////////////#

## libfile ##
target = libutils.a

# complier and options
FC = ifort
FLAG = -O3 -fpp -fPIC
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

# directories
MyLib = ../

# compile links
MODD = -I$(MyLib)/mod
OPTION = $(MODD)

#//// files to be compiled ////#
# independent modules
OBJ1 = const.o funcs.o random.o io.o utils.o array.o 
OBJ2 = readfile.o 
OBJ3 = myfunc.o cmblbisp.o stat.o mycls.o mycmb.o mygal.o myGW.o mySNeIa.o myrec.o 
OBJ = ${OBJ1} ${OBJ2} ${OBJ3}

$(target): $(OBJ)
	ar crs $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	mv *.mod ${MyLib}/mod; mv $(target) ${MyLib}/lib/

clean:
	rm -f *.o* *.mod *.a

# file dependencies

myfunc.o : funcs.o utils.o
cmblbisp.o : const.o utils.o myfunc.o
stat.o : array.o utils.o
mycls.o : utils.o const.o
mycmb.o : readfile.o
mygal.o : utils.o funcs.o
myGW.o  : const.o myfunc.o
myrec.o : readfile.o

