! The Wolff algorithm for the XY model.
! A clock version of the XY model is used. 
! Program written by Jos Thijssen, 2007
! Program discussed in chapter 15 of `Computational Physics',
! second edition, Cambridge 2007

PROGRAM XY

  IMPLICIT NONE
  INTEGER, PARAMETER :: Size=60, MaxClock = 64
  INTEGER :: Lattice(Size, Size), Step, MaxStep
  DOUBLE PRECISION :: JCoupl, CosArr(MaxClock), SinArr(MaxClock), &  
                      CosArr2(2*MaxClock), PI, HelMod, &
                      TotCos, TotSinX, TotSinY
  LOGICAL Tagged(Size, Size)

  CALL Initialise
  MaxStep = 30000
! Equilibrate
  DO Step = 1, 100
    CALL WolffStep
    IF (MODULO(Step, 30) == 0) CALL Display()
  END DO
! Harvest results
  DO Step = 1, MaxStep
    CALL WolffStep
    CALL CalcHelMod
    IF (MODULO(Step, 100) == 0) CALL Display()
  END DO
! Calculate helicity modulus
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
! Fill arays with the necessary values of sin and cos
  DO I=1, MaxClock
    CosArr(I) = COS(2*PI*(I-1)/DBLE(MaxClock))
    CosArr2(2*I-1) = COS(2*PI*(2*I-2)/DBLE(2*MaxClock))
    CosArr2(2*I) = COS(2*PI*(2*I-1)/DBLE(2*MaxClock))
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
      X2 = X1 + 0.75D0*COS(Lattice(IX,IY)*2.D0*PI/MaxClock)
      Y2 = Y1 + 0.75D0*SIN(Lattice(IX,IY)*2.D0*PI/MaxClock)
      CALL Draw(X1, Y1, X2, Y2)
    END DO
  END DO
#endif
  END SUBROUTINE Display

  SUBROUTINE WolffStep
! The Wolff algorithm.
  INTEGER :: I, J, U
  DOUBLE PRECISION :: RI, RJ, RClock
  
  Tagged = .FALSE.
  CALL Random_Number(RI)
  I = INT(RI*Size)+1
  CALL Random_Number(RJ)
  J = INT(RJ*Size)+1
  CALL Random_Number(RClock)
! Note that U should assume integer or half-integer values. 
! If only integer values would be allowed, then an initial
! configuration with all spins 0 would lead to only even spins.
! Our U is is an integer twice as big as its `actual' value, which therefore 
! can be half-integer
  U = INT(2*RClock*MaxClock)
  CALL GrowCluster(I, J, U)
  END SUBROUTINE WolffStep


  RECURSIVE SUBROUTINE GrowCluster(I, J, U)
  INTEGER, INTENT(IN) :: I, J, U
  INTEGER :: IP1, IM1, JP1, JM1, PrevSpin

  Tagged(I,J) = .TRUE.
  ! Flip the spin. Note that the integer U is twice the value of 
  ! the angle corresponding to the vector used in the reflection step
  Lattice(I,J) = MODULO(U-Lattice(I,J)+MaxClock/2, MaxClock)
  IP1=MODULO(I, Size)+1
  IM1=MODULO(I-2, Size)+1
  JP1=MODULO(J, Size)+1
  JM1=MODULO(J-2, Size)+1
  PrevSpin = Lattice(I,J)
  IF (.NOT.(Tagged(IM1, J))) CALL TryAdd(IM1, J, PrevSpin, U)
  IF (.NOT.(Tagged(IP1, J))) CALL TryAdd(IP1, J, PrevSpin, U)
  IF (.NOT.(Tagged(I, JM1))) CALL TryAdd(I, JM1, PrevSpin, U)
  IF (.NOT.(Tagged(I, JP1))) CALL TryAdd(I, JP1, PrevSpin, U)
  END SUBROUTINE GrowCluster




  RECURSIVE SUBROUTINE TryAdd(I, J, PrevSpin, U)
  INTEGER, INTENT(IN) :: I, J, PrevSpin, U
  DOUBLE PRECISION :: R, Pf, X

  ! Calculate the new J_ij. Note that the `finer' array CosArr2 is 
  ! used instead of CosArr as the actual value U/2 may be half-integer.
  X = 2*JCoupl*CosArr2(MODULO(2*PrevSpin-U,2*MaxClock)+1)*&
               CosArr2(MODULO(2*Lattice(I,J)-U, 2*MaxClock)+1)
!  Oops! The book says EXP(-X) rather than X. Sorry!
  Pf = 1.D0-EXP(X)
  CALL Random_Number(R)
  IF (R<Pf) THEN
    CALL GrowCluster(I, J, U)
  END IF
  END SUBROUTINE TryAdd






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
