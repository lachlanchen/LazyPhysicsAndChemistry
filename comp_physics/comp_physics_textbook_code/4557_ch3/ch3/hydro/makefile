# COMPHY directory path:
COMPHY	=	$(HOME)/comphy
LDFLAGS	=	-O
FCOPTS	=	-O -c 
FC	=	gfortran
LIBS	=	$(COMPHY)/lib/geneig.a -llapack

%.o:%.f
	$(FC) $(FCOPTS) $<

Hatom: Hatom.o globHatom
	$(FC) $(LDFLAGS) -o Hatom Hatom.o $(LIBS)
clean:
	rm -f Hatom.o Hatom
