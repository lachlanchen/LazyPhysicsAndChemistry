# Makefile to compile Quantum Monte Carlo simulation for
# Computational Physics

# COMPHY directory path:
COMPHY	=	$(HOME)/comphy
LIBS	=	-L$(COMPHY)/lib -lran
FFLAGS	=	-O
LDFLAGS	=	-O
FC	=	gfortran

ho: ho.o
	$(FC) -o ho $(LDFLAGS) ho.o $(LIBS)

ho.o: ho.f globHOFP
	$(FC) $(FFLAGS) -c ho.f
clean:
	rm -f *.o ho
