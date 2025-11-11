      IMPLICIT NONE

      INTEGER MaxSize, MaxHist
      PARAMETER (MaxSize = 250, MaxHist=50)

      INTEGER Size, MaxStep, Histogram(-MaxHist:MaxHist), HistSize,
     .        InitStep

      REAL*8 Chain(0:MaxSize,3), Epsilon, Energy, InvEps, 
     .       MaxDisp, TotE, PotArray(0:49,0:49,0:20)

      COMMON Chain, Epsilon, MaxDisp, InvEps, Energy, TotE,
     .       PotArray, Size, MaxStep, Histogram,
     .       HistSize, InitStep
