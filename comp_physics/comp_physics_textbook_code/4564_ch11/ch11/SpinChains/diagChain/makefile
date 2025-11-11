FC 	= 	ifort
LD	=	ifort
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


heisen: heisen.o lanczos.o
	$(FC) $(LDOPTS) heisen.o lanczos.o -o heisen $(LIBS)

symheis: symheis.o
	$(FC) $(LDOPTS) symheis.o -o symheis $(LIBS)

clean:
	rm -f *.mod *.o  *~ heisen symheis
