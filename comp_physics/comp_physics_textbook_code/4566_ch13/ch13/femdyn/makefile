FC 	= 	gfortran
LD	=	gfortran
FCOPTS	=	-O -c 
LDOPTS	=   	-O

LIBS=-L/usr/X11R6/lib -L/home/jost/XPS/lib -llapack -lxps -lX11
%.mod: %.f90
	$(FC) $(FCOPTS) $<
%.o : %.f90
	$(FC) $(FCOPTS) $<
%.o : %.f
	$(FC) -c $<
%.o : %.c
	cc -c -I/usr/X11R6/include/X11/ $<


femdyn: femdyn.o
	$(FC) $(LDOPTS) femdyn.o -o femdyn $(LIBS)
clean:
	rm -f *.o *.mod *~ femdyn
