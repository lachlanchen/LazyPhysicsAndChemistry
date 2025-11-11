FC 	= 	gfortran
LD	=	gfortran
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


sw: SW.o
	$(FC) $(LDOPTS) SW.o -o sw $(LIBS)
bt: BackTrack.o
	$(FC) $(LDOPTS) BackTrack.o -o bt $(LIBS)
wolff: Wolff.o
	$(FC) $(LDOPTS) Wolff.o -o wolff $(LIBS)
wolff_xy: wolff_xy.o
	$(FC) $(LDOPTS) wolff_xy.o -o wolff_xy $(LIBS)
xy: xy.o
	$(FC) $(LDOPTS) xy.o -o xy $(LIBS)
clean:
	rm -f *.o *.mod *~ bt sw wolff wolff_xy xy
