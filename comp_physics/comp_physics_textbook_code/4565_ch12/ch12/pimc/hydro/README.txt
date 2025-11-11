The program hydro needs a file CoulPot which contains 
the cumulant potential. This potential depends on the "time step"
Delta tau. The potential can be generated using the program
"gen_cumulant" which is found in the subdir "cumulant".
This program must be run with the same value of Delta tau
as the program hydro.
You need many steps for obtaining good statistics. A
sample input file is provided. Just run the program as
./hydro < input


