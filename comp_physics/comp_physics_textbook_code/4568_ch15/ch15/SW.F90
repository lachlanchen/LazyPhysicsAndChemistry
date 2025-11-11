! Program for the Swendsen-Wang algorithm. The Hoshen-Kopelman 
! algorithm is used to identify the different clusters.
! Program written by Jos Thijssen, 2007

PROGRAM Swendsen_Wang

  IMPLICIT NONE
  INTEGER, PARAMETER :: Size=60, MaxLabel = 10000, MaxClust = 3000
  INTEGER :: Lattice (Size, Size), Links(Size,Size,2), Roots(Size, Size), &
             Step, MaxStep
  INTEGER :: LL(MaxLabel), ClustCnt, ClusterIndex(MaxClust), ClusterSpin(MaxClust)
  DOUBLE PRECISION :: FreezeFac, Cv, JCoupl, E, TotE, TotESq

  CALL Initialise
  MaxStep = 6000
  TotE = 0
  TotESq = 0
! Equilibration
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
! Find clusters using the Hishen Kopelman algorithm
  CALL IdentifyClusters
! Adapt the values at each site...
  CALL SetRoots
! Now assign the (random) ClusterSpin values to each cluster.
  DO I=1, Size
    DO J=1, Size
      R = Roots(I,J)
      Lattice(I,J) = ClusterSpin(R)
    END DO
  END DO
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
! Hoshen-Kopelman
  INTEGER :: I, J, IP1, JP1, IM1, JM1, NeighCnt, Help, Help2, &
             KK, Props(4), TempProp
  LOGICAL :: NewProp

! LL is a sort of linked list. Each element of this list points 
! to the next as long as LL(R) < 0. When LL(R)>0, R is the 'proper label'
! and LL(R) is the cluster size (not really necessary for SW)
  LL = 0
  Roots = 0
  ClustCnt = 0
  DO I=1, Size
    IP1=MODULO(I,Size)+1
    IM1=MODULO(I-2,Size)+1
    DO J=1, Size
!      print *, 'l', I, J, LL(2)
      JP1=MODULO(J,Size)+1
      JM1=MODULO(J-2,Size)+1
      NeighCnt = 0
      ! Check if the first neighbour has already been assigned to a cluster
      IF ((Roots(IM1,J)/=0) .AND. (Links(IM1,J,1)/=0)) THEN
        IF (Lattice(IM1, J)*Lattice(I,J)<0) print *, 'a hallo', I, J
        NeighCnt = NeighCnt + 1
        CALL Classify(IM1, J, Props(NeighCnt))
      END IF
      ! Check whether the 2nd neighbour has already been assigned to a cluster
      ! and whether this is different from the previous (if any)
      IF ((Roots(IP1,J)/=0) .AND. &
          (Links(I,J,1)/=0) .AND. (I==Size)) THEN
        CALL Classify(IP1, J, TempProp)
        NewProp = .TRUE.
        DO KK=1, NeighCnt
          IF (TempProp == Props(KK)) NewProp = .FALSE.
        END DO
        IF (NewProp) THEN
          NeighCnt = NeighCnt + 1
          Props(NeighCnt) = TempProp
        END IF
      END IF
      ! Check whether the 3rd neighbour has already been assigned to a cluster
      ! and whether this is different from the previous (if any)
      IF ((Roots(I,JM1)/=0) .AND. (Links(I,JM1,2)/=0)) THEN
        CALL Classify(I, JM1, TempProp)
        NewProp = .TRUE.
        DO KK=1, NeighCnt
          IF (TempProp == Props(KK)) NewProp = .FALSE.
        END DO
        IF (NewProp) THEN
          NeighCnt = NeighCnt + 1
          Props(NeighCnt) = TempProp
        END IF
      END IF
      ! Check whether the 4th neighbour has already been assigned to a cluster
      ! and whether this is different from the previous (if any)
      IF ((Roots(I,JP1)/=0) .AND. &
          (Links(I,J,2)/=0) .AND. (J==Size)) THEN
        CALL Classify(I, JP1, TempProp)
        NewProp = .TRUE.
        DO KK=1, NeighCnt
          IF (TempProp == Props(KK)) NewProp = .FALSE.
        END DO
        IF (NewProp) THEN
          NeighCnt = NeighCnt + 1
          Props(NeighCnt) = TempProp
        END IF
      END IF
      ! A new cluster is assigned to this site
      IF (NeighCnt == 0) THEN
        ClustCnt = ClustCnt + 1
        Roots(I,J) = ClustCnt
        LL(ClustCnt) = 1
      ELSE
        ! The lowest propval, 'Help', of the neighbours is identified.
        ! The different cluster sizes of the neighbours are added to the
        ! current one (1) and LL(Help) points to this value (total cluster size)
        ! The roots(I,J) at the current site is given the value Help.
        ! The neighbouring propvals all point to Help
        Help = MinVal(Props(1:NeighCnt))
        IF (Help<0) print *, 'problem: help < 0!!'
        Roots(I,J) = Help
        Help2 = 0
        DO KK=1, NeighCnt
!          print *, i, j, LL(Props(KK))
          Help2 = Help2+LL(Props(KK))
          LL(Props(KK)) = -Help
        END DO
        LL(Help) = Help2 + 1
      END IF
    END DO
  END DO
  END SUBROUTINE IdentifyClusters


  SUBROUTINE Classify(I, J, PropVal)
! This routine returns the proper value of the 
! cluster site i, j is part of
  INTEGER :: I, J, PropVal, T, R

  R = Roots(I,J)
  T = -LL(R)
  IF (T<0) THEN
    PropVal = R
  ELSE
    R = T
    T = -LL(R)
    IF (T<0) THEN
      PropVal = R
      LL(Roots(I,J)) = -PropVal
    ELSE
      DO WHILE(T>=0)
        R = T
        T = -LL(R)
      END DO
      PropVal = R
      LL(Roots(I,J)) = -PropVal
    END IF
  END IF
  END SUBROUTINE Classify




  SUBROUTINE SetRoots
! Adjust the roots to assume the proper value of 
! cluster the sites are part of.
  INTEGER :: I, J, T, R, Val
  DOUBLE PRECISION :: Rand

  ClustCnt = 0
  ClusterIndex = 0
  ClusterSpin = 1
  DO I=1, Size
    DO J=1, Size
      R = Roots(I,J)
      IF (R /= 0) THEN
        T = -LL(R)
        DO WHILE(T>0)
          R = T
          T = -LL(R)
        END DO
        Roots(I,J) = R
      ELSE
        print *, 'root is nul!!!', I, J
      END IF
      IF (ClusterIndex(R) == 0) THEN
        ClustCnt = ClustCnt+1
        ClusterIndex(R) = ClustCnt
!        print *, i, j, clustercount
        CALL Random_Number(Rand)
        IF (Rand<0.5D0) ClusterSpin(R) = -1
      END IF
    END DO
  END DO
  END SUBROUTINE SetRoots

END PROGRAM Swendsen_Wang