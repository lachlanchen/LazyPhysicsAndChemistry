! Program for simulating the flow of an incompressible fluid using
! Lattice Boltzmann simulation. The fluid runs through a 2D pipe and
! is blocked by a rectangular object.
! Program described in `Computational Physics', second edition, 
! J. M. Thijssen, Cambridge 2007
! Program written 2003-2007



PROGRAM LB
 IMPLICIT NONE

 INTEGER, PARAMETER :: Lx=120, Ly=70, DIM=2, VelNum=9
 DOUBLE PRECISION :: Densities(Lx, Ly, VelNum), NewDens(Lx, Ly, VelNum), &
                     AverVel(Lx,Ly,DIM), AverDens(Lx, Ly), TauInv
 INTEGER :: Alive(Lx, Ly)
 INTEGER :: BasVec(VelNum,DIM)

 CALL Initialise
 CALL Simulate
 CALL EndProg

CONTAINS

SUBROUTINE Initialise
  INTEGER :: I, J, K
  DOUBLE PRECISION :: InitDens, EqDist(VelNum)
! Alive denotes whether a site is blocked or not
  Alive = 0
  Alive(:,1) = 4
  Alive(:,Ly) = 4
  Alive(40:80,17:51) = 4
  InitDens = 1.0D0
  AverVel = 0.D0
  AverDens = 1.D0
  DO I=1, Lx
    DO J=1, Ly
      CALL CalcEquiDist(I, J, EqDist)
      DO K=1, VelNum
        Densities(I,J,K) = EqDist(K)
        NewDens(I,J,K) = Densities(I,J,K)
      END DO
    END DO
  END DO
  NewDens = Densities
! TauInv is the inverse relaxation time, which determines the 
! viscusity nu (see book)
  TauInv = 1.25D0
  CALL FillBasVec
  CALL CalcAverVel

  CALL StartGraphics
END SUBROUTINE Initialise


SUBROUTINE FillBasVec
! Orientations on d2q9 lattice
  INTEGER :: K
  BasVec(1,1) = 0
  BasVec(1,2) = 0
  DO K=1, 4
    BasVec(2*K,1) = MOD(2-K,2)
    BasVec(2*K,2) = MOD(3-K,2)
  END DO
  BasVec(3,1) = 1
  BasVec(3,2) = 1
  BasVec(5,1) = -1
  BasVec(5,2) = 1
  BasVec(7,1) = -1
  BasVec(7,2) = -1
  BasVec(9,1) = 1
  BasVec(9,2) = -1
  DO K=1, VelNum
    print *, basvec(k,1), basvec(k,2)
  end do
END SUBROUTINE FillBasVec



SUBROUTINE StartGraphics
#ifdef Plot
  CALL InitPlot("lightblue", 1200, 300, "out.ps", 1)
  CALL Framing (-1.0D0, -1.0D0, Lx+1.0D0, Ly+1.0d0)
  CALL PutStopButton
  CALL PutContButton
  CALL Draw(0.D0, 0.D0, DBLE(Ly), 0.D0)
#endif
END SUBROUTINE StartGraphics




SUBROUTINE EndProg
! Writing and plotting final results
  INTEGER :: Num, R, G, B, X, Y, J, I
  DOUBLE PRECISION :: U
  OPEN (8, File='veloc.dat')
  OPEN (9, FILE='dens.dat')
! Various velocity and density profiles are written to files
  write (8,'(I8,2F12.8)') ((J, AverVel(Lx/2, J, 1), AverVel(Lx/2, J, 2)), J=1, Ly)
  write (9,*) ((J, AverDens(Lx/2, J)), J=1, Ly)
  write (8,*)
  write (9,*)
  write (8,'(I8,2F12.8)') ((J, AverVel(Lx/3, J, 1), AverVel(Lx/3, J, 2)), J=1, Ly)
  write (9,*) ((J, AverDens(Lx/3, J)), J=1, Ly)
  write (8,*)
  write (9,*)
  write (8,'(I8,2F12.8)') ((J, AverVel(Lx/10, J, 1), AverVel(Lx/10, J, 2)), J=1, Ly)
  write (9,*) ((J, AverDens(Lx/10, J)), J=1, Ly)
  write (8,*)
  write (9,*)
  write (8,'(I8,2F12.8)') ((I, AverVel(I, Ly/2, 1), AverVel(I, Ly/2, 2)),I=1, Lx)
  write (9,*) ((I, AverDens(I, Ly/2)), I=1, Lx)
  CLOSE(8)
  CLOSE(9)
#ifdef Plot
  CALL InitPlot('white', Lx*10,Ly*10,'pressure.ps',1)
  CALL Framing (0.D0, 0.D0, 1.D0, 1.D0)
  CALL PutStopButton()
  DO y = 1, Ly 
    DO x = 1, Lx 
      u = sqrt(dot_product(AverVel(x,y,:),AverVel(x,y,:)))
      Num = MIN(INT(U*1000*255),255)
      r = Num; g = int (Num*Num/256.0); b = 255-Num
      IF (Alive(x,y) .NE.0) THEN
        r=0
        g=0
        b=0
      END IF
      CALL setcol(r, g, b)
      CALL SetPSColor(dble(r/256.d0), dble(g/256.d0), &
                      dble(b/256.d0))
      CALL FillRectangle(dble(x)/dble(Lx), dble(y)/dble(Ly),&
                        dble(x+1)/dble(Lx), dble(y+1)/dble(Ly))
    END DO
  END DO

  CALL EndPlot()
#endif
END SUBROUTINE EndProg



SUBROUTINE Simulate
! Main simulation routine
  INTEGER :: I, J, K, Time
  INTEGER, PARAMETER :: EndTime=2000

  CALL CalcAverVel
  DO Time = 1, EndTime
    IF (MOD(Time, 100).EQ.0) THEN
      print *, 'Step', Time
      CALL DrawConfig
    END IF
  CALL FillNewDensities
  END DO
END SUBROUTINE Simulate


SUBROUTINE DrawConfig
  INTEGER :: I, J, K
  DOUBLE PRECISION :: Vx, Vy, X1, X2, Y1, Y2
#ifdef Plot
  CALL SetNamedBackground('lightblue');
  DO I=1, Lx, 1
    DO J=1, Ly, 1
      IF (Alive(I,J).NE. 4) THEN
        vx = SUM(BasVec(:,1)*Densities(I,J,:))
        vy = SUM(BasVec(:,2)*Densities(I,J,:))
        X1 = DBLE(I)
        Y1 = DBLE(J)
        X2 = X1+200*Vx
        Y2 = Y1+200*Vy
        CALL Draw(X1, Y1, X2, Y2)
      END IF
    END DO
  END DO
#endif
END SUBROUTINE DrawConfig


SUBROUTINE CalcAverVel
! Calculates the average velocity on each site, which is then 
! used as input for calculating the equilibrium velocity at that
! velocity
  INTEGER I, J, K
  DOUBLE PRECISION X1, Y1, X2, Y2, DeltaP, DeltaV, IP, DeltaIP, &
                   UDU, TotDens

  DeltaV = 0.002D0
  TotDens = 0.D0
  DO I=1, Lx
    DO J=1, Ly
      AverVel(I,J,1) = SUM(BasVec(:,1)*Densities(I,J,:))
      AverVel(I,J,2) = SUM(BasVec(:,2)*Densities(I,J,:))
    END DO
  END DO
  TotDens = SUM(AverDens)
  AverVel(:,:,1) = AverVel(:,:,1)/AverDens
  AverVel(:,:,1) = AverVel(:,:,1)/AverDens
END SUBROUTINE CalcAverVel



SUBROUTINE CalcEquiDist(I,J, EqDist)
! Calculate equilibrium distribution using the given average velocity

  INTEGER I, J, K
  DOUBLE PRECISION EqDist(VelNum), SqVel, F0, IP, &
                   X1, X2, Y1, Y2, Vx, Vy

  SqVel = AverVel(I,J,1)**2+AverVel(I,J,2)**2
  EqDist(1) = 4*AverDens(I,J)/9.D0*(1-1.5d0*SqVel)
  DO K=2, VelNum
    IP = AverVel(I,J,1)*BasVec(K,1)+AverVel(I,J,2)*BasVec(K,2)
    IF (Mod(K,2).EQ.0) THEN
      F0 = AverDens(I,J)/9.D0
    ELSE
      F0 = AverDens(I,J)/36.D0
    END IF
    EqDist(K) = F0*(1.D0+3*IP-1.5D0*SqVel+4.5D0*IP*IP)
  END DO
END SUBROUTINE CalcEquiDist






SUBROUTINE FillNewDensities
! Replace the old densities by the (partly relaxed) new ones.
! This routine also calls BoundaryConditionsNew which takes into
! account the BC at the boundaries of the pipe and the blocking 
! object. 
! Also, DeltaV is added to the horizontal velocity

  INTEGER I, J, K, NewI, NewJ, NewK
  DOUBLE PRECISION EqDist(VelNum), DeltaV

  DeltaV = 0.000001D0
  CALL BoundaryConditionsNew
  CALL CalcAverVel
  DO I=1, Lx
    DO J=1, Ly
      IF (Alive(I,J).EQ.0) THEN
        CALL CalcEquiDist(I, J, EqDist)
        NewDens(I,J,:)=Densities(I,J,:)*(1-TauInv)+EqDist(:)*TauInv
      END IF
    END DO
  END DO
  WHERE (Alive==0)
    NewDens(:,:,2) = NewDens(:,:,2)+DeltaV
    NewDens(:,:,6) = NewDens(:,:,6)-DeltaV
  END WHERE
END SUBROUTINE FillNewDensities





SUBROUTINE BoundaryConditionsNew

  INTEGER I, J, K, NewI, NewJ, NewK
  DOUBLE PRECISION LocDens

  DO I=1, Lx
    DO J=1, Ly
      IF (Alive(I, J) .NE. 0) THEN
        Densities(I,J,1) = NewDens(I, J, 1)
        DO K=2, VelNum
          NewI = MOD(I-1+BasVec(K,1)+Lx,Lx)+1
          NewJ = J+BasVec(K,2)
          NewK = MOD(K-2+4, 8) + 2
          IF ((NewJ<=Ly .AND. NewJ>=1).AND.&
                  (Alive(NewI, NewJ).NE.-4)) THEN
                Densities(I,J,K) = NewDens(NewI, NewJ, NewK)
          END IF
        END DO
      END IF
    END DO
  END DO
  DO I=1, Lx
    DO J=1, Ly
      IF (Alive(I,J).NE.0) THEN
        NewDens(I,J,:) = Densities(I,J,:)
      END IF
    END DO
  END DO
  Densities(:,2:Ly-1,:) = 0.D0
  DO I=1, Lx
    DO J=1, Ly
      Densities(I,J,1) = NewDens(I,J,1)
      DO K=2, VelNum
        NewI = MOD(I-1+BasVec(K,1)+Lx,Lx)+1
        NewJ = J+BasVec(K,2)
        IF (NewJ.GT.1 .AND. NewJ.LT.Ly) THEN
          NewK = K
          Densities(NewI,NewJ,NewK) = NewDens(I,J,K)
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE BoundaryConditionsNew



END PROGRAM LB