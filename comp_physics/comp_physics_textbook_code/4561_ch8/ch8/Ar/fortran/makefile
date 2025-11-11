# COMPHY directory path:
COMPHY	=	$(HOME)/comphy
FC 	= 	gfortran
FCOPTS	=	-O -c
LDOPTS	=	-O
LIBS  	=	-L$(COMPHY)/lib -lran

md: md.o
	$(FC) $(LDOPTS) -o md md.o $(LIBS)
md.o: md.f md.glb
	$(FC) $(FCOPTS) md.f
clean:
	rm -f *.o md
