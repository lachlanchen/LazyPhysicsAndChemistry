FC 	= 	gfortran
LD	=	gfortran
FCOPTS	=	-O -c 
LDOPTS	= 	-O

LIBS=-L/usr/X11R6/lib -L/home/jost/XPS/lib -llapack -lxps -lX11
%.mod: %.f90
	$(FC) $(FCOPTS) $<
%.o : %.f90
	$(FC) $(FCOPTS) $<
%.o : %.f
	$(FC) -c $<
%.o : %.c
	cc -c -I/usr/X11R6/include/X11/ $<

rivara: rivardat.o rivara.o triangfem.o mainrivar.o
	$(FC) $(LDOPTS) rivardat.o rivara.o triangfem.o mainrivar.o -o rivara $(LIBS) 

clean:
	rm -f *.o *.mod *~ rivara
