FC 	= 	ifort 
LD	=	ifort 
FCOPTS	=	-O -c -DPlot
LDOPTS	= -O
COMPHY  = $(HOME)/comphy
XPSLIBS = -L/usr/X11R6/lib -L$(XPSDIR) -lxps -lX11
# Uses the XPS library. XPSDIR is an environment variable containing the 
# path to the "lib" subdir of the XPS directory.
# If you do not use graphics, change -DPlot into -UPlot above,
# and remove $(XPSLIBS) from the line below.
LIBS= -L$(COMPHY)/lib $(XPSLIBS) -lran
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


md: md.o
	$(LD) $(LDOPTS) md.o -o md $(LIBS)

clean:
	rm -f *.mod *.o  *~ md
