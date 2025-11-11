FC 	= 	gfortran
LD	=	gfortran
FCOPTS	=	-c -O
LDOPTS	= 	-O

LIBS=-L/usr/X11R6/lib -L/home/jost/XPS/lib \
 -llapack -lxps -lX11
%.mod: %.f90
	$(FC) $(FCOPTS) $<
%.o : %.f90
	$(FC) $(FCOPTS) $<
%.o : %.f
	$(FC) -c $<
%.o : %.c
	cc -c -I/usr/X11R6/include/X11/ $<

main: brmod.o md2d.o main.o rangen.o
	$(FC) $(LDOPTS) main.o brmod.o md2d.o rangen.o -o main $(LIBS) 
brmod.o: brmod.f90
	$(FC) $(FCOPTS) brmod.f90 
md2d.o: md2d.f90
	$(FC) $(FCOPTS) md2d.f90  
main.o: main.f90
	$(FC) $(FCOPTS) main.f90 

clean: 
	rm -rf *.o main 

