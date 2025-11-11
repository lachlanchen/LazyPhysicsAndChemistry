# Makefile to compile Quantum Monte Carlo simulation for
# Computational Physics

# COMPHY directory path:
COMPHY	=	$(HOME)/comphy
LIBS	=	-L$(COMPHY)/lib -lran
# Add -lm to LIBS on HP-UX 10
FFLAGS	=	-O -c
LDFLAGS	=	-O
FC	=	gfortran
# fort77 for HP-UX 10

hodmc: hodmc.o
	$(FC) -o hodmc $(LDFLAGS) hodmc.o $(LIBS)

hodmc.o: hodmc.f globHODMC
	$(FC) $(FFLAGS) hodmc.f
clean:
	rm -f *.o hodmc
