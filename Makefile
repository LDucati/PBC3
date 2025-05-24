#

### Fortran compiler
FC     = gfortran
FC90   = gfortran

LINKER = gfortran

### Disable/Enable OpenMP support
OMP = -fopenmp

### Compilation/Linking options
FOPT = -g -Wall -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow 

FTUNE = -O3 -ffast-math 

### Extra libraries
#FLIB=-llapack
FLIB =  

FFLAGS = $(OMP)

.SUFFIXES: .f .F .F90 .f90

BIN = pbc3.x

OBJ = 

OBJ90 = geompar.o 


.F.o:
	$(FC) $(FOPT) $(FTUNE) $(OMP) -c $< -o $@
.F90.o: 
	$(FC90) $(FOPT) $(FTUNE) $(OMP) -c $< -o $@
.f90.o: 
	$(FC90) $(FOPT) $(FTUNE) $(OMP) -c $< -o $@

all: $(BIN)

pbc3.x: $(OBJ) $(OBJ90) main.o
	$(LINKER) -o pbc3.x main.o $(OBJ) $(OBJ90) $(FFLAGS) $(FLIB)

clean:
	rm -f $(BIN) *.o *.mod



