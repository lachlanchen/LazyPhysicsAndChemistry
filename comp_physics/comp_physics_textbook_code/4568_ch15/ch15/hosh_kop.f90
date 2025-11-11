! Straightforward implementation and check of the
! Hoshen-Kopelman algorithm

PROGRAM Hosh_Kop

  IMPLICIT NONE
  INTEGER, PARAMETER :: Size=10, MaxLabel = 100
  INTEGER :: Lattice (Size, Size), Links(Size,Size,2), Roots(Size, Size)
  INTEGER :: LL(MaxLabel), ClustCnt

  CALL InitLattice
  CALL IdentifyClusters
  CALL PrintLat
  print '(10I5)', Roots

CONTAINS

  SUBROUTINE InitLattice
  INTEGER :: I, J, IP1, JP1, IM1, JM1
  DOUBLE PRECISION :: R

  Lattice = -1
  Links = 0
  DO I=1, Size
    DO J=1, Size
      CALL Random_Number(R)
      IF (R<0.5D0) Lattice(I,J) = 1
    END DO
  END DO
  DO I=1, Size
    DO J=1,Size
      IP1=MODULO(I,Size)+1
      JP1=MODULO(J,Size)+1
      IF (Lattice(I,J) == Lattice(IP1, J)) THEN 
        CALL Random_Number(R)
        IF (R<1.d5) Links(I, J, 1) = 1
      END IF
      IF (Lattice(I,J) == Lattice(I, JP1)) THEN 
        CALL Random_Number(R)
        IF (R<1.d5) Links(I, J, 2) = 1
      END IF
    END DO
  END DO
  Roots = 0
  ClustCnt = 0
  END SUBROUTINE InitLattice


  SUBROUTINE IdentifyClusters
  INTEGER :: I, J, IP1, JP1, IM1, JM1, NeighCnt, Help, Help2, KK, Props(4)

  DO I=1, Size
    DO J=1, Size
      IF (Lattice(I,J) == 1) THEN
        IP1=MODULO(I,Size)+1
        JP1=MODULO(J,Size)+1
        IM1=MODULO(I-2,Size)+1
        JM1=MODULO(J-2,Size)+1
        NeighCnt = 0
        IF ((Lattice(IM1,J)==1) .AND. (Roots(IM1,J)/=0) .AND. (Links(IM1,J,1)/=0)) THEN
          NeighCnt = NeighCnt + 1
          CALL Classify(IM1, J, Props(NeighCnt))
        END IF
        IF ((Lattice(IP1,J)==1) .AND. (Roots(IP1,J)/=0) .AND. &
            (Links(I,J,1)/=0) .AND. (I==Size)) THEN
          NeighCnt = NeighCnt + 1
          CALL Classify(IP1, J, Props(NeighCnt))
        END IF

        IF ((Lattice(I,JM1)==1) .AND. (Roots(I,JM1)/=0) .AND. (Links(I,JM1,2)/=0)) THEN
          NeighCnt = NeighCnt + 1
          CALL Classify(I, JM1, Props(NeighCnt))
        END IF
        IF ((Lattice(I,JP1)==1) .AND. (Roots(I,JP1)/=0) .AND. &
            (Links(I,J,2)/=0) .AND. (J==Size)) THEN
          NeighCnt = NeighCnt + 1
          CALL Classify(I, JP1, Props(NeighCnt))
        END IF
        IF (NeighCnt == 0) THEN
          ClustCnt = ClustCnt + 1
          Roots(I,J) = ClustCnt
          LL(ClustCnt) = 1
          print *, i, j, roots(i,j), NeighCnt
        ELSE
          Help = MinVal(Props(1:NeighCnt))
          Roots(I,J) = Help
          Help2 = 0
          DO KK=1, NeighCnt
            Help2 = Help2+LL(Props(KK))
            LL(Props(KK)) = -Help
          END DO
          LL(Help) = Help2 + 1
        END IF
      END IF
    END DO
  END DO
  END SUBROUTINE IdentifyClusters


  SUBROUTINE Classify(I, J, PropVal)
  INTEGER :: I, J, PropVal, T, R

  R = Roots(I,J)
  T = -LL(R)
  IF (T<0) THEN
    PropVal = R
    return
  END IF
  DO WHILE(T>0)
    print *, 'a', t
    R = T
    T = -LL(R)
  END DO
  PropVal = R
  LL(Roots(I,J)) = -PropVal
  END SUBROUTINE Classify

  SUBROUTINE PrintLat
  INTEGER :: I, J, T, R, Val

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
      END IF
    END DO
  END DO
  END SUBROUTINE PrintLat

END PROGRAM Hosh_Kop