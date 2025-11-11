FC = 		gfortran
FCOPTS = 	-O -c
LDOPTS =	-O 
%.o : %.f
	$(FC) $(FCOPTS) $<

gencumulant: gencumulant.o simpson.o
	$(FC) $(LDOPTS) -o gencumulant gencumulant.o simpson.o

clean:
	rm -f *.o gencumulant