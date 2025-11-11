      IMPLICIT NONE

      INTEGER MaxSize, MaxHist
      PARAMETER (MaxSize = 250, MaxHist=75)

      INTEGER Size, MaxStep, Histogram(-MaxHist:MaxHist), 
     .        InitStep, HistSize

      REAL*8 Chain(0:MaxSize), Epsilon, Energy, InvEps, 
     .       MaxDisp, TotE

      COMMON Chain, Epsilon, MaxDisp, InvEps, Energy, TotE,
     .       Size, MaxStep, Histogram, 
     .       HistSize, InitStep
