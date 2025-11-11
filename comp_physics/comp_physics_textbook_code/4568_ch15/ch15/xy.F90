! Plain Metropolis MC for the XY model
! Program written by Jos Thijssen, 2007
! See chapter 15 of the book `Computational Physics'
! by Jos Thijssen. Cambridge, second edition, 2007

PROGRAM XY

  IMPLICIT NONE
  INTEGER, PARAMETER :: Size=30, MaxClock = 64
  INTEGER :: Lattice(Size, Size), Step, MaxStep
  DOUBLE PRECISION :: CosArr(MaxClock), SinArr(MaxClock), &  
                      PI, HelMod, TotCos, TotSinX, TotSinY, JCoupl
  LOGICAL Tagged(Size, Size)

  CALL Initialise
  MaxStep = 6000 ! Number of production steps
  DO Step = 1, 100
    CALL MCStep
    IF (MODULO(Step, 30) == 0) CALL Display()
  END DO
  DO Step = 1, MaxStep
    CALL MCStep
    CALL CalcHelMod
    IF (MODULO(Step, 100) == 0) CALL Display()
  END DO
  TotCos  = TotCos/MaxStep
  TotSinX = TotSinX/MaxStep
  TotSinY = TotSinY/MaxStep
  HelMod  = (TotCos - JCoupl*(TotSinX+TotSinY))/(2*Size*Size)
  print *, HelMod, TotCos, TotSinX
#ifdef Plot
  CALL EndPlot()
#endif
CONTAINS

  SUBROUTINE Initialise
  INTEGER :: I

  PI = 4.D0*ATAN(1.D0)
  PRINT *, 'Give coupling constant J/kT'
  READ *, JCoupl
  DO I=1, MaxClock
    CosArr(I) = COS(2*PI*(I-1)/DBLE(MaxClock))
    SinArr(I) = SIN(2*PI*(I-1)/DBLE(MaxClock))
  END DO
#ifdef Plot
  CALL InitPlot('lightblue', 700,700, 'out.ps', 1)
  CALL Framing(-0.05D0*Size, -0.05D0*Size, &
                1.05D0*Size,  1.05D0*Size)
  CALL PutStopButton()
#endif
  Lattice = 0
  TotCos = 0.D0
  TotSinX = 0.D0
  TotSinY = 0.D0
  END SUBROUTINE Initialise



  SUBROUTINE Display

  INTEGER :: IX, IY
  DOUBLE PRECISION :: X1, Y1, X2, Y2
#ifdef Plot
  CALL NextPage()
  DO IX=1, Size
    DO IY=1, Size
      X1 = DBLE(IX)
      Y1 = DBLE(IY)
      X2 = X1 + 0.6D0*COS(Lattice(IX,IY)*2.D0*PI/MaxClock)
      Y2 = Y1 + 0.6D0*SIN(Lattice(IX,IY)*2.D0*PI/MaxClock)
      CALL Draw(X1, Y1, X2, Y2)
    END DO
  END DO
#endif
  END SUBROUTINE Display



  SUBROUTINE MCStep
  INTEGER, PARAMETER :: MaxDev = MaxClock/6
  INTEGER :: I, J, IP1, IM1, JP1, JM1, OldSpin, NewSpin, Dev
  DOUBLE PRECISION :: OldEner, NewEner, BoltzFac, R

  DO I=1, Size
    IP1 = MODULO (I, Size) + 1
    IM1 = MODULO (I-2, Size) + 1
    DO J=1, Size
      JP1 = MODULO (J, Size) + 1
      JM1 = MODULO (J-2, Size) + 1
      OldSpin = Lattice(I,J);
      CALL Random_Number(R);
      Dev = INT(2*R*MaxDev+1)-MaxDev
      NewSpin = OldSpin + Dev
      OldEner = CosArr(MODULO(OldSpin-Lattice(IM1,J), MaxClock)+1)
      NewEner = CosArr(MODULO(NewSpin-Lattice(IM1,J), MaxClock)+1)
      OldEner = OldEner+CosArr(MODULO(OldSpin-Lattice(IP1,J), MaxClock)+1)
      NewEner = NewEner+CosArr(MODULO(NewSpin-Lattice(IP1,J), MaxClock)+1)
      OldEner = OldEner+CosArr(MODULO(OldSpin-Lattice(I,JM1), MaxClock)+1)
      NewEner = NewEner+CosArr(MODULO(NewSpin-Lattice(I,JM1), MaxClock)+1)
      OldEner = OldEner+CosArr(MODULO(OldSpin-Lattice(I,JP1), MaxClock)+1)
      NewEner = NewEner+CosArr(MODULO(NewSpin-Lattice(I,JP1), MaxClock)+1)
      BoltzFac = Exp(JCoupl*(NewEner-OldEner))
      CALL Random_Number(R);
      IF (R<BoltzFac) Lattice(I,J) = NewSpin
    END DO
  END DO
  END SUBROUTINE MCStep




  SUBROUTINE CalcHelMod
  INTEGER :: I, J, IP1, JP1, DiffAng, SpinVal
  DOUBLE PRECISION :: LocTotSinX, LocTotSinY
! See Eq. (15.97). TotCos and TotSin are used above in the calculation
! of HelMod

  LocTotSinX = 0.D0
  LocTotSinY = 0.D0
  DO I=1, Size
    IP1 = MODULO (I, Size) + 1
    DO J=1, Size
      JP1 = MODULO (J, Size) + 1
      SpinVal = Lattice(I, J)
      DiffAng = MODULO(SpinVal -Lattice(IP1, J), MaxClock) + 1
      TotCos = TotCos + CosArr(DiffAng)
      LocTotSinX = LocTotSinX + SinArr(DiffAng)
      DiffAng = MODULO(SpinVal -Lattice(I, JP1), MaxClock) + 1
      TotCos = TotCos + CosArr(DiffAng)
      LocTotSinY = LocTotSinY + SinArr(DiffAng)
    END DO
  END DO
  TotSinX = TotSinX + LocTotSinX**2
  TotSinY = TotSinY + LocTotSinY**2
  END SUBROUTINE CalcHelMod

END PROGRAM XY
