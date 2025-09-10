#
# This program may need lapack and blas routines
#
# For MPI support, we would also need MPI libs and include files
#
# Note: we define two compiler variables: $FC for fixed format F77 or F90 and $FC90 for
# free format F90


### Fortran compiler
FC     = gfortran
FC90   = gfortran
LINKER = gfortran

### Disable/Enable OpenMP support
OMP = 

### Compilation/Linking options
FOPT = -O2 -g -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow # gfortran 

### Extra libraries
#FLIB = -llapack
FLIB =

FFLAGS = $(OMP)

.SUFFIXES: .f .F .F90 .f90

BIN = parse-namelists

OBJ = 

OBJ90 = definitions.o  \
        constants-parameters.o \
	molecule-basis-module.o \
	namelists.o \
        factorials.o \
        check-print-basis.o \
	overlap.o \
        bas-powers.o


.F.o:
	$(FC) $(FOPT) -c $< -o $@
.F90.o: 
	$(FC90) $(FOPT) -c $< -o $@
.f90.o: 
	$(FC90) $(FOPT) -c $< -o $@

all: $(BIN)

parse-namelists: $(OBJ) $(OBJ90) parse-namelists.o
	$(LINKER) -o parse-namelists parse-namelists.o $(OBJ) $(OBJ90) $(FFLAGS) $(FLIB)

clean:
	rm -f $(BIN) *.o *.mod



