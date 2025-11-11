# This programs works in its current form only for the ifort compiler.
# Use gf_fem.f90 for the gfortran compiler
 

FC 	= 	ifort
LD	=	ifort
FCOPTS	=	-O -c 
LDOPTS	= 	-O
XPSLIBS = -L$(XPSDIR)/lib -lxps -L/usr/X11R6/lib -lxps -lX11

LIBS= $(XPSLIBS) -llapack 
%.mod: %.f90
	$(FC) $(FCOPTS) $<
%.o : %.f90
	$(FC) $(FCOPTS) $<
%.o : %.f
	$(FC) -c $<
%.o : %.c
	cc -c -I/usr/X11R6/include/X11/ $<


fem: fem.o
	$(FC) fem.o -o fem $(LIBS)
gf_fem: gf_fem.o
	$(FC) gf_fem.o -o gf_fem $(LIBS)

clean:
	rm -f *.mod *.o  *~ fem gf_fem
