! Program for the Swendsen-Wang algorithm. The back-track 
! algorithm is used to identify the different clusters.
! Program written by Jos Thijssen, 2007

PROGRAM Swendsen_Wang

  IMPLICIT NONE
  INTEGER, PARAMETER :: Size=60, MaxClust = 3000
  INTEGER :: Lattice (Size, Size), Links(Size,Size,2), &
             TotMag, Step, MaxStep
  !INTEGER :: ClustCnt, ClusterIndex(MaxClust), ClusterSpin(MaxClust)
  DOUBLE PRECISION :: FreezeFac, Cv, JCoupl, E, TotE, TotESq
  LOGICAL Tagged(Size, Size)

  CALL Initialise
  MaxStep = 6000
  TotE = 0
  TotESq = 0
  ! Equilibrate
  DO Step = 1, 100
    CALL SWStep
  END DO
  ! Production
  DO Step = 1, MaxStep
    CALL SWStep
    CALL CalcEnergies
    IF (MODULO(Step, 100) == 0) CALL Display()
  END DO
  Cv = JCoupl*JCoupl*(TotESq/DBLE(MaxStep)-(TotE/DBLE(MaxStep))**2)/(Size*Size)
  print *, cv
#ifdef Plot
  CALL EndPlot()
#endif
CONTAINS

  SUBROUTINE Initialise

  PRINT *, 'Give coupling constant J/kT'
  READ *, JCoupl
  FreezeFac= EXP(-2*JCoupl)
#ifdef Plot
  CALL InitPlot('lightblue', 700,700, 'out.ps', 1)
  CALL Framing(-0.05D0*Size, -0.05D0*Size, &
                1.05D0*Size,  1.05D0*Size)
  CALL PutStopButton()
#endif
  Lattice = 1
  END SUBROUTINE Initialise



  SUBROUTINE Display

  INTEGER :: IX, IY
#ifdef Plot
  CALL NextPage()
  DO IX=1, Size
    DO IY=1, Size
      IF (Lattice(IX,IY) .EQ. 1) THEN
        CALL FillRectangle(DBLE(IX), DBLE(IY), &
                           DBLE(IX+1), DBLE(IY+1))
!        print *, IX, IY
      END IF
    END DO
  END DO
#endif
  END SUBROUTINE Display




  SUBROUTINE SWStep
! Straightforward SW algorithm. The difficulty is in IdentifyClusters!
  INTEGER :: I, J, IP1, JP1
  DOUBLE PRECISION :: R
  
  Links = 0
  DO I=1, Size
    IP1=MODULO(I,Size)+1
    DO J=1, Size
      JP1=MODULO(J,Size)+1
! Freeze some of the links according to the SW algorithm
      IF (Lattice(I,J) == Lattice(IP1, J)) THEN
        CALL Random_Number(R)
        IF (R>FreezeFac) Links(I,J,1) = 1
      END IF
      IF (Lattice(I,J) == Lattice(I, JP1)) THEN
        CALL Random_Number(R)
        IF (R>FreezeFac) Links(I,J,2) = 1
      END IF
    END DO
  END DO
! Find clusters using the BackTrack algorithm
  CALL IdentifyClusters
  END SUBROUTINE SWStep




  SUBROUTINE CalcEnergies
  INTEGER :: I, J, IP1, JP1

  E = 0
  DO I=1, Size
    IP1=MODULO(I,Size)+1
    DO J=1, Size
      JP1=MODULO(J,Size)+1
      E = E - Lattice(I, J)*(Lattice(IP1, J) + Lattice(I, JP1))
    END DO
  END DO
  TotE = TotE + E
  TotESq = TotESq + E**2
  END SUBROUTINE CalcEnergies



  SUBROUTINE IdentifyClusters

  INTEGER :: I, J
  DOUBLE PRECISION :: R
  
  Tagged = .FALSE.
  DO I=1, Size
    DO J=1, Size
      IF (.NOT. Tagged(I, J)) THEN
        CALL Random_Number(R)
        IF (R<0.5D0) THEN
          CALL BackTrack(I,J, 1)
        ELSE
          CALL BackTrack(I,J, -1)
        END IF
      END IF
    END DO
  END DO
  END SUBROUTINE IdentifyClusters

  RECURSIVE SUBROUTINE BackTrack (I, J, ClusterSpin)
  INTEGER, INTENT(IN) :: I, J, ClusterSpin
  INTEGER :: IP1, IM1, JP1, JM1

  Tagged(I,J) = .TRUE.
  Lattice(I,J) = ClusterSpin
  IP1=MODULO(I, Size)+1
  IM1=MODULO(I-2, Size)+1
  JP1=MODULO(J, Size)+1
  JM1=MODULO(J-2, Size)+1
  IF (.NOT.(Tagged(IM1, J)).AND.(Links(IM1,J,   1)/=0)) CALL BackTrack(IM1, J, ClusterSpin)
  IF (.NOT.(Tagged(IP1, J)).AND.(Links(I,  J,   1)/=0)) CALL BackTrack(IP1, J, ClusterSpin)
  IF (.NOT.(Tagged(I, JM1)).AND.(Links(I,  JM1, 2)/=0)) CALL BackTrack(I, JM1, ClusterSpin)
  IF (.NOT.(Tagged(I, JP1)).AND.(Links(I,  J,   2)/=0)) CALL BackTrack(I, JP1, ClusterSpin)
  END SUBROUTINE BackTrack
  
END PROGRAM Swendsen_Wang
