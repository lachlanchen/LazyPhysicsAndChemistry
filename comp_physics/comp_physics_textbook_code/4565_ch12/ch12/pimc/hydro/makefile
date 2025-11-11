FC 	= 	gfortran
LD	=	gfortran
FCOPTS	=	-g -c -DPlot
LDOPTS	= 	-g
COMPHY  =	$(HOME)/comphy
XPSLIBS = 	-L$(COMPHY)/lib -lran -L$(XPSDIR)/lib -lxps -L/usr/X11R6/lib -lxps -lX11
LIBS= $(XPSLIBS) -llapack 
%.o : %.f
	$(FC) $(FCOPTS) $<
%.o : %.F
	$(FC) $(FCOPTS) $<


hydro: hydro.o
	$(FC) hydro.o -o hydro $(LIBS)

clean:
	rm -f *.mod *.o  *~ hydro
