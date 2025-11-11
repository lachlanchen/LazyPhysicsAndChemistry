      PROGRAM CalcRes
C This program reads in the results for the Virial 
C (in reduced units of epsilon), and the temperature 
C (in units of epsilon/k_Boltzmann) and the potential
C (in units of epsilon) and the density (in 1/sigma^3)
C The progrm then returns the values for
C potential (including the tail beyond r_cut-off) and
C the quantity P/(k_Boltzmann T).

      IMPLICIT NONE

      DOUBLE PRECISION Temp, Virial, Pot, RMin,
     .                 PkT, TotalPot, PI, Density
      INTEGER N

      PI = 4.D0*DATAN(1.D0)

      PRINT *, 'Give Temperature, Virial Potential,'
      PRINT *, 'Density, N and R_cut-off'
      READ *, Temp, Virial, Pot, Density, N, RMin
      TotalPot = Pot/N - 8*PI*Density/(3*RMin**3)
      PkT = 1.D0-Virial/(3*N*Temp)-8*PI*Density/(3*Temp*2.5**3)
      PRINT *, 'Total Potential = ', TotalPot
      PRINT *, 'P/kT = ', PkT
      END
