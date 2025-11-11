# COMPHY directory path:
COMPHY	=	$(HOME)/comphy
XPSDIR	=	$(HOME)/XPS/lib
FC 	= 	gfortran
XPSLIBS	=	-L/usr/X11R6/lib -L$(XPSDIR) -lxps -lX11 	
FCOPTS	=	-O -DPlot -c
# The Plot flag enables plotting using the XPS package.
# This package can also be found in the computational
# physics program directory
# To switch this option off, use -UPlot, and remove $(XPSLIBS) below
LDOPTS	=	-O
LIBS	=	-L$(COMPHY)/lib  -lran  $(XPSLIBS)

nitro: nitro.o
	$(FC) $(LDOPTS) -o nitro nitro.o $(LIBS)
nitro.o: nitro.F nitro.glob
	$(FC) $(FCOPTS) nitro.F

clean:
	rm -rf *.o nitro
