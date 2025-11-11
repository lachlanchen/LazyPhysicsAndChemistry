# Makefile to compile Quantum Monte Carlo simulation for
# Computational Physics

# COMPHY directory path:
COMPHY	=	$(HOME)/comphy
LIBS	=	-L$(COMPHY)/lib -lran
FFLAGS	=	-O -c
LDFLAGS	=	-O
FC	=	gfortran
# fort77 for HP-UX 10

ho: ho.o
	$(FC) -o ho $(LDFLAGS) ho.o $(LIBS)

ho.o: ho.f globHO
	$(FC) $(FFLAGS) ho.f
clean:
	rm -f *.o ho
