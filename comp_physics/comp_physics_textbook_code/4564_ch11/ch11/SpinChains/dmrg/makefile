FC 	= 	gfortran
LD	=	gfortran
FCOPTS	=	-O -c 
LDOPTS	=	-O

LIBS= -llapack
%.mod: %.f90
	$(FC) $(FCOPTS) $<
%.o : %.f90
	$(FC) $(FCOPTS) $<
%.o : %.f
	$(FC) -c $<
%.o : %.c
	cc -c $<


white: white.o
	$(FC) $(LDOPTS) white.o -o white $(LIBS)

noack: noack.o
	$(FC) $(LDOPTS) noack.o -o noack $(LIBS)

clean:
	rm -f *.mod *.o  *~ white noack
