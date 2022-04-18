FC = ifort
FLAG = -qopenmp -ip -fpp
FLAGS = $(FLAG) $(DBAG)
OPT =  -I${HOME}/Work/Lib/mylib/mod -I${HOME}/Work/Lib/mylib/pub/Healpix/include  -L${HOME}/Work/Lib/mylib/lib -L${HOME}/Work/Lib/mylib/pub/Healpix/lib -L${HOME}/Work/Lib/mylib/pub/cfitsio  -lnldd -lanafull -lutils -lhealpix -lcfitsio
OBJ = main.o
all: exe
%.o: %.f90
	 $(FC) $(FLAGS) $(OPT) -c $*.f90
exe: $(OBJ)
	 $(FC) $(FLAGS) $(OBJ) $(OPT) -o $@
clean:
	 rm -f *.o* *.e* *.mod *.d *.pc *.obj core* *.il
