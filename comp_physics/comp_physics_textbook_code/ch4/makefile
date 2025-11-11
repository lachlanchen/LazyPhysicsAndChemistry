# COMPHY directory path:
COMPHY	=	$(HOME)/comphy
LDFLAGS	=	-O  
FC	=	gfortran
LIBS	=	$(COMPHY)/lib/geneig.a -llapack
FFLAGS 	= 	-O -c

%.o:%.f
	$(FC) $(FFLAGS) $<

he: he.o globHe
	$(FC) $(LDFLAGS) -o he he.o $(LIBS)

h2: h2.o globH2
	$(FC) $(LDFLAGS) -o h2 h2.o $(LIBS)

opthe: opthe.o globHe
	$(FC) $(LDFLAGS) -o opthe opthe.o $(LIBS)

clean:
	rm -f *.o opthe he h2
