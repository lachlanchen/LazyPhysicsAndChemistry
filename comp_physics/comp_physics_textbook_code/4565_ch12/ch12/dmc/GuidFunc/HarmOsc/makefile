# Makefile to compile Quantum Monte Carlo simulation for
# Computational Physics

# COMPHY directory path:
COMPHY	=	$(HOME)/comphy
LIBS	=	-L$(COMPHY)/lib -lran
FFLAGS	=	-O -c
LDFLAGS	=	-O
FC	=	gfortran

hodmc: hodmc.o
	$(FC) -o hodmc $(LDFLAGS) hodmc.o $(LIBS)

hodmc.o: hodmc.f globHODMC
	$(FC) $(FFLAGS) -c hodmc.f
clean:
	rm -f *.o hodmc
