# COMPHY directory path:
COMPHY	=	$(HOME)/comphy
LDFLAGS	=	-O
FCOPTS	=	-O -c 
FC	=	gfortran
LIBS	=	$(COMPHY)/lib/geneig.a -llapack

%.o:%.f
	$(FC) $(FCOPTS) $<

deepwell: deepwell.o globvarc
	$(FC) $(LDFLAGS) -o deepwell deepwell.o $(LIBS)

deepeven: deepeven.o globvarc
	$(FC) $(LDFLAGS) -o deepeven deepeven.o $(LIBS)

deepodd: deepodd.o globvarc
	$(FC) $(LDFLAGS) -o deepodd deepodd.o $(LIBS) 

clean:
	rm -f *.o deepwell deepodd deepeven
