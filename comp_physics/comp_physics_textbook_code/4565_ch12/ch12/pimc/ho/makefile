# This programs works in its current form only for the ifort compiler.
# Use gf_fem.f90 for the gfortran compiler
 

FC 	= 	ifort
LD	=	ifort
FCOPTS	=	-g -c 
LDOPTS	= 	-g
COMPHY	=	/home/jost/comphy
# XPSLIBS = -L$(XPSDIR)/lib -lxps -L/usr/X11R6/lib -lxps -lX11

LIBS= -L$(COMPHY)/lib -lran
%.mod: %.f90
	$(FC) $(FCOPTS) $<
%.o : %.f90
	$(FC) $(FCOPTS) $<
%.o : %.f
	$(FC) $(FCOPTS) -c $<
%.o : %.c
	cc -c -I/usr/X11R6/include/X11/ $<


ho: ho.o
	$(FC) $(LDOPTS) ho.o -o ho $(LIBS)
ho3: ho3.o
	$(FC) $(LDOPTS) ho3.o -o ho3 $(LIBS)

clean:
	rm -f *.mod *.o  *~ ho
