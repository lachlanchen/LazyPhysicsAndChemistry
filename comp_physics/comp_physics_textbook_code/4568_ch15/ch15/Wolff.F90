! Program for the Wolff algorithm. Ising version. 
! Program written by Jos Thijssen, 2007
! Program described in `Computational Physics'
! by Jos Thijssen, second edition, Cambridge 2007

PROGRAM Wolff

  IMPLICIT NONE
  INTEGER, PARAMETER :: Size=60
  INTEGER :: Lattice (Size, Size), Links(Size,Size,2), &
             Step, MaxStep
! Cv is the specific heat, which is calculated from Eav and Esq_av
  DOUBLE PRECISION :: FreezeFac, Cv, JCoupl, E, TotE, TotESq
  LOGICAL Tagged(Size, Size)

  CALL Initialise
  MaxStep = 6000
  TotE = 0
  TotESq = 0
  DO Step = 1, 100
! Equilibration
    CALL WolffStep
    IF (MODULO(Step, 30) == 0) CALL Display()
  END DO
  DO Step = 1, MaxStep
    CALL WolffStep
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




  SUBROUTINE WolffStep
! A call is made to GrowCluster, starting from a random 
! position
  INTEGER :: I, J, ClusterSpin
  DOUBLE PRECISION :: RI, RJ
  
  Tagged = .FALSE.
  CALL Random_Number(RI)
  I = INT(RI*Size)+1
  CALL Random_Number(RJ)
  J = INT(RJ*Size)+1
  ClusterSpin = - Lattice(I,J)
  CALL GrowCluster(I, J, ClusterSpin)
  END SUBROUTINE WolffStep

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


  RECURSIVE SUBROUTINE GrowCluster(I, J, ClusterSpin)
! Cluster is grown by trying to add neighbours of present sites
  INTEGER, INTENT(IN) :: I, J, ClusterSpin
  INTEGER :: IP1, IM1, JP1, JM1

  Tagged(I,J) = .TRUE.
  Lattice(I,J) = ClusterSpin
  IP1=MODULO(I, Size)+1
  IM1=MODULO(I-2, Size)+1
  JP1=MODULO(J, Size)+1
  JM1=MODULO(J-2, Size)+1
  IF (.NOT.(Tagged(IM1, J))) CALL TryAdd(IM1, J, ClusterSpin)
  IF (.NOT.(Tagged(IP1, J))) CALL TryAdd(IP1, J, ClusterSpin)
  IF (.NOT.(Tagged(I, JM1))) CALL TryAdd(I, JM1, ClusterSpin)
  IF (.NOT.(Tagged(I, JP1))) CALL TryAdd(I, JP1, ClusterSpin)
  END SUBROUTINE GrowCluster

  RECURSIVE SUBROUTINE TryAdd(I, J, ClusterSpin)
! We try to add the spin at I,J to the Wolff cluster
  INTEGER, INTENT(IN) :: I, J, ClusterSpin
  DOUBLE PRECISION :: R

  IF (Lattice(I,J)*ClusterSpin == -1) THEN
    CALL Random_Number(R)
    IF (R>FreezeFac) THEN
      CALL GrowCluster(I, J, ClusterSpin)
    END IF
  END IF
  END SUBROUTINE TryAdd



  
END PROGRAM Wolff