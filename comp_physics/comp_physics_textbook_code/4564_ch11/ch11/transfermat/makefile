COMPHY	=	$(HOME)/comphy
LDFLAGS	=	-O -C
LIBS	=       $(COMPHY)/lib/libran.a -llapack
FC = gfortran
FFLAGS = -O -C

tm: tm.o
	$(FC) $(LDFLAGS) tm.o -o tm $(LIBS)


clean:
	rm -rf *.o tm