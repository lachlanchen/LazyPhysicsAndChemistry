# Makefile to compile Quantum Monte Carlo simulation for
# Computational Physics

# COMPHY directory path:
COMPHY	=	$(HOME)/comphy
LIBS	=	-L$(COMPHY)/lib -lran
FFLAGS	=	-O -c
LDFLAGS	=	-O
FC	=	gfortran

hydr: hydr.o
	$(FC) -o hydr $(LDFLAGS) hydr.o $(LIBS)

hydr.o: hydr.f globHydrDMC
	$(FC) $(FFLAGS) -c hydr.f
clean:
	rm -f *.o hydr
