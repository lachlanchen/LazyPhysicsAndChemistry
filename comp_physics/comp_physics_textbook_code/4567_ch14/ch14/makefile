FC 	= 	ifort
LD	=	ifort
# Change -DPlot into -UPlot if you do not use graphics
FCOPTS	=	-O -c -DPlot
LDOPTS	= 	-O

# Remove line below if you do not run graphics package
LIBS=-L/usr/X11R6/lib -L/home/jost/XPS/lib -lxps -lX11
%.mod: %.f90
	$(FC) $(FCOPTS) $<
%.o : %.f90
	$(FC) $(FCOPTS) $<
%.o : %.F90
	$(FC) $(FCOPTS) $<
%.o : %.f
	$(FC) -c $<
%.o : %.c
	cc -c -I/usr/X11R6/include/X11/ $<


lb: lb.o setcol.o
	$(FC) lb.o setcol.o -o lb $(LIBS)
clean:
	rm -f *.o *.mod *~ lb
