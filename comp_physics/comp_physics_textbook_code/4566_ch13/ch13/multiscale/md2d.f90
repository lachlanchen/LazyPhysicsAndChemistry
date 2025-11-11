MODULE MD2D
!*****************************************************************
!* This program performs a Molecular Dynamics (MD) simulation    *
!* for a monatomic material.                                     *
!* Periodic boundary conditions are assumed. The potential is    *
!* cut off beyond a distance "Rcutoff". Furthermore, the minimum *
!* image convention is used to calculate forces etc.             *
!*                                                               *
!*                                                               *
!* Program written by Jos Thijssen	                         *
!* Program written 2004                                          *
!*****************************************************************
  IMPLICIT NONE

  INTEGER, PARAMETER :: Dim=2

  DOUBLE PRECISION :: &
            Temperature, Volume, InitTime, SimTime, TimeStep, &
            VolSize(DIM), CorrStep, ROuter, RInner, Sep, Fac2

  DOUBLE PRECISION, ALLOCATABLE :: Pos(:,:), Momentum(:,:), PairDist(:,:), MDForce(:,:), &
                                   FEMPart(:,:), FEMForce(:,:)


  INTEGER :: ScaleTime, PairNum, DispInt, Nx, Ny, PartNum, MaxLen, &
             InitStep, PairListTime

  INTEGER, ALLOCATABLE :: PairList(:,:), BorderPart(:)

  LOGICAL :: Scale


CONTAINS




SUBROUTINE InitMD
  CALL InitParameters
  CALL InitPositions
  CALL InitMomenta
END SUBROUTINE InitMD





SUBROUTINE InitParameters

  INTEGER :: I, Num

  OPEN (UNIT=8, FILE='md2d.in')

  READ (8,*) Temperature
  READ (8,*) TimeStep
  READ (8,*) InitTime
  READ (8,*) ScaleTime
  READ (8,*) SimTime
  READ (8,*) VolSize
  READ (8,*) Ny
  READ (8,*) num
  READ (8,*) DispInt

  InitStep = NINT(InitTime/TimeStep)
  PairListTime = 10000 ! Large enough to have very few changes of the pairlist
!****** Density and particle number define the volume of the system **
  CorrStep = 0.02D0
  CALL InitRand(num)
  RInner = 2.5D0
  ROuter = 3.3D0
  Nx = NINT(4*Ny*VolSize(1)/VolSize(2))
  PartNum = Nx*Ny
  ALLOCATE (Pos(PartNum, DIM))
  ALLOCATE (MDForce(PartNum, DIM))
  ALLOCATE (Momentum(PartNum, DIM))
  MaxLen = PartNum**2/4
  ALLOCATE (PairDist(MaxLen,4))
  ALLOCATE (PairList(MaxLen,2))
  ALLOCATE (FEMPart(2*Ny, DIM))
  ALLOCATE (FEMForce(2*Ny, DIM))
  ALLOCATE (BorderPart(Ny))
  VolSize(1) = VolSize(1)*5 ! We want to have extra empty space left and right
END SUBROUTINE InitParameters



SUBROUTINE InitPositions

! *** Positions are stored on a regular fcc lattice. Lattice constant* 
! *** is adjusted to volume size to fill the volume homogeneously ****

INTEGER :: IX, IY, IZ, Counter
DOUBLE PRECISION :: LattConst=0.996891D0, YLattConst, XOffset, ScaleConst

!*******Calculate Volume Size from Volume ************************
!*******LinCell is number of cells along one side ****************
  WRITE (6,*) VolSize, PartNum
  ScaleConst = 2**(1.D0/6)
  LattConst = LattConst*ScaleConst
  YLattConst = LattConst*0.5D0*SQRT(3.D0)
  print *, lattconst, ylattconst, ny, 'piep'

  XOffset = 24.d0*0.996891D0*ScaleConst!1.12246D0
  Counter = 0
  print *, 'in md2, Ny is ', Ny, PartNum, Nx*Ny
  DO IX = 0, Nx-1
    DO IY = 0, Ny-1
      Counter = Counter + 1
      IF (MOD(IY,2).EQ.0) THEN
        Pos(Counter,:) = (/IX*LattConst+XOffset, IY*YLattConst/)
      ELSE
        Pos(Counter,:) = (/(IX+0.5D0)*LattConst+XOffset, IY*YLattConst/)
      END IF
      IF (IX == 0) THEN
        BorderPart(IY+1) = Counter
      END IF
    END DO
    print *, 'check', SUM(Pos(Counter-Ny+1:Counter,2)-(Ny-1)*YLattConst*0.5)/Ny
  END DO
  print *, 'borderpart', borderpart
  DO IY=1, Ny
    print *, IY, Pos(borderpart(IY),:)
  END DO
END SUBROUTINE InitPositions



SUBROUTINE InitMomenta

!     Initialise momenta of particles. All velocity components
!     are drawn from random generator with Gaussian distribution

  INTEGER :: I

  DOUBLE PRECISION :: R1, R2, TotEner, V(DIM)

! *****Assign initial velocities to all particles*****************
! *****TotEner is used for rescaling the velocities***************
  TotEner = 0
! Specifically for 3D!
  DO I = 1, PartNum
    Call ExpRand(R1, R2)
    Momentum(I,:) = (/R1, R2/)
  END DO

! *****Set Centre of Mass velocity equal to zero*******************
  V = 0
  DO I=1,PartNum
    V = V + Momentum(I,:)
  END DO
  V = V/PartNum

  DO I=1,PartNum
    Momentum(I,:) = Momentum(I,:) - V
  END DO

! ***** Rescale velocities to the right temperature ***************
  Call CalcTemp (TotEner)
  Call Rescale (TotEner)
END SUBROUTINE InitMomenta



SUBROUTINE InitSim (DeltaT)

! *******************************************************************
! ****  Subroutine in which the actual simulation is performed. *****
! **** IF "Scale" is .TRUE., velocities are regularly rescaled ******
! **** "Time" defines simulation time. ******************************
! **** Forces of previous integration step are stored in OldForceX **
! **** etc. ********************************************************* 
! *******************************************************************

  DOUBLE PRECISION ::  &
      TotEner, Virial, DeltaT

  INTEGER :: StepNum, Step
  LOGICAL :: Scale, Special, VirCalc
  Scale = .FALSE.
  TimeStep = DeltaT 
  CALL UpdatePairList
  CALL CalcPairList
  CALL CalcMDForce (MDForce)

!  StepNum = INT(Time/TimeStep)
  Fac2 = 0.5D0*TimeStep*TimeStep
  Momentum = 0.D0
  CALL CalcMDForce(MDForce)
END SUBROUTINE InitSim

SUBROUTINE DoMDStep1(Step, DeltaT)
  DOUBLE PRECISION :: Aver1, Aver2, MaxVal, ScaleConst, &
             DeltaT, Virial, TotEner, LattConst, YLattConst, XOffset
  LOGICAL :: Scale, Special, VirCalc
  INTEGER :: Step, I, K, IX, IY, Ctr, Ct1, Ct2, MaxPos
  TimeStep = DeltaT
  Pos = Pos + TimeStep*Momentum + Fac2*MDForce
  IF (MOD(Step, DispInt)==0) THEN
!    CALL SetNamedBackground('lightblue')
!    CALL SetNamedColor('red')
!    DO I=1, PartNum
!      CALL SetPoint(Pos(I,1), Pos(I,2))
!    END DO
    ScaleConst = 2**(1.D0/6)
    LattConst = ScaleConst*0.996891D0
    YLattConst = LattConst*0.5*SQRT(3.D0)
    XOffset = 24.D0*LattConst
    CALL SetNamedColor('black')
    CALL SetPoint(XOffSet-0.5D0, 0.D0)
    Ctr = 0
    MaxPos = 0
    MaxVal = 0.D0
    DO IX = 0, Nx-1
      Aver1 = 0.D0
      Aver2 = 0.D0
      Ct1=0; Ct2=0
      DO IY = 0, Ny-1
        Ctr = Ctr + 1
        IF (MOD(IY,2).EQ.0) THEN
          Aver1 = Aver1 + Pos(Ctr,2)-IY*YLattConst
          Ct1 = Ct1 + 1
        ELSE
          Ct2 = Ct2 + 1
          Aver2 = Aver2 + Pos(Ctr,2)-IY*YLattConst
        END IF
      END DO
      IF ((Aver1/Ct1 >MaxVal).OR.(Aver2/Ct2>MaxVal)) THEN 
        MaxVal = MAX(Aver1/Ct1, Aver2/Ct2)
        MaxPos = IX
      END IF
!      IF (IX==0) THEN
!        CALL SetPoint(DBLE(IX)*1.12246D0+XOffset-0.5D0, 10000*(Aver1+Aver2)/(Ct1+Ct2))
!      ELSE
        CALL DrawTo(DBLE(IX)*1.12246D0+XOffset-0.5D0, 10000*(Aver1+Aver2)/(Ct1+Ct2))
!      END IF
    END DO
    IF (MOD(Step, 1000)==0) THEN
      CALL CalcElast()
    END IF
  END IF
  Momentum = Momentum + 0.5D0*TimeStep*MDForce
END SUBROUTINE DoMDStep1

SUBROUTINE DoMDStep2(Step, DeltaT)
  INTEGER :: Step
  DOUBLE PRECISION :: DeltaT, Virial, TotEner

  Momentum = Momentum + 0.5D0*TimeStep*MDForce
  IF (MOD(Step, PairListTime)==0) THEN
    CALL UpdatePairList
  END IF
END SUBROUTINE DoMDStep2





SUBROUTINE CopyFEMToMD(FemPos, N)
  INTEGER :: N
  DOUBLE PRECISION :: FemPos(N,Dim)

  FEMPart(:,:) = FemPos(:,:)
END SUBROUTINE CopyFEMToMD

SUBROUTINE OutputMDBorder(MDPos, N)
  INTEGER :: N, I
  DOUBLE PRECISION :: MDPos(N,Dim)

  DO I=1, Ny
    MDPos(I,:) = Pos(BorderPart(I),:)
  END DO
END SUBROUTINE OutputMDBorder



SUBROUTINE CalcMDForce(MDForce)

! *** Calculation of forces. We assume that the forces are a super- **
! *** position of central-symmetric forces between two particles.   ** 

  DOUBLE PRECISION :: MDForce(PartNum,Dim),&
                    R2, R4, R8, R14, ForceConst, S(DIM), &
                    RMin2, Virial, aver

  INTEGER :: I, J, k, PairCnt

  MDForce = 0.D0
  Virial = 0.D0
  CALL CalcPairList
  DO PairCnt = 1,PairNum
    I = PairList(PairCnt,1)
    J = PairList(PairCnt,2)
    R2 = PairDist(PairCnt,4)
    IF (R2.LT.Rinner*RInner) THEN
      S = PairDist(PairCnt,1:DIM)
      RMin2 = 1/R2
      R4 = RMin2*RMin2
      R8 = R4*R4
      R14= R8*R4*RMin2

      ForceConst = 24.d0*(2*R14-R8)
      MDForce(I,:) = MDForce(I,:) + S*ForceConst

      MDForce(J,:) = MDForce(J,:) - S*ForceConst
    END IF
  END DO
  FEMForce = 0.D0
  DO I = 1, Ny
    DO J=1, 2*Ny
      S =  Pos(BorderPart(I),:)-FEMPart(J,:)
      R2 = DOT_PRODUCT(S,S)
      IF (R2.LT.Rinner*RInner) THEN
        RMin2 = 1/R2
        R4 = RMin2*RMin2
        R8 = R4*R4
        R14= R8*R4*RMin2
        ForceConst = 24.d0*(2*R14-R8)
        MDForce(BorderPart(I),:) = MDForce(BorderPart(I),:) + 0.5D0*S*ForceConst
        FEMForce(J,:) = FEMForce(J,:) - 0.5D0*S*ForceConst
      END IF
    END DO
  END DO
END SUBROUTINE CalcMDForce


SUBROUTINE CalcElast
! *** Calculation of the elastic matrix **

  DOUBLE PRECISION :: C(Dim, Dim, Dim, Dim),&
                    R2, R4, R8, R10, R14, R16, ElastConst, S(DIM), &
                    RMin2, Virial, C3(3,3), akappa, alambda, anu, amu, &
                    ElastC2

  INTEGER :: I, J, PairCnt, kappa, lambda, mu, nu

  C = 0.D0
  DO PairCnt = 1,PairNum
    I = PairList(PairCnt,1)
    J = PairList(PairCnt,2)
    R2 = PairDist(PairCnt,4)
    IF (R2.LT.Rinner*RInner) THEN
      S = PairDist(PairCnt,1:DIM)
      RMin2 = 1/R2
      R4 = RMin2*RMin2
      R8 = R4*R4
      R10 = R8*RMin2
      R14= R8*R4*RMin2
      R16 = R8*R8
    
      ElastConst = (672*R16-192*R10)
      ElastC2 = -48*r14+24*r8
      DO kappa = 1, DIM
        akappa = s(kappa)
        DO lambda = 1, DIM
          alambda = s(lambda)
          DO mu = 1, DIM
            amu = s(mu)
            DO nu = 1, DIM
              anu = s(nu)
              C(kappa, lambda, mu, nu) = C(kappa, lambda, mu, nu) + &
                                    ElastConst*akappa*alambda*amu*anu
              IF (Lambda == nu) THEN
                C(kappa, lambda, mu, nu) = C(kappa, lambda, mu, nu) + &
                    ElastC2*akappa*amu
              END IF
            END DO
          END DO
        END DO
      END DO
    END IF
  END DO
  C = C/(PartNum*0.5*sqrt(3.D0)*(1.12246D0*0.99681D0)**2)
  C3(1,1) = C(1,1,1,1); C3(2,2) = C(2,2,2,2)
  C3(1,2) = C(1,1,2,2); C3(2,1) = C(2,2,1,1)
  C3(1,3) = 0.5D0*(C(1,1,1,2)+C(1,1,2,1)); 
  C3(3,1) = 0.5D0*(C(1,2,1,1)+C(2,1,1,1));
  C3(2,3) = 0.5D0*(C(2,2,1,2)+C(2,2,2,1)); 
  C3(3,2) = 0.5D0*(C(1,2,2,2)+C(2,1,2,2));
  C3(3,3) = 0.25D0*(C(1,2,1,2)+C(2,1,2,1)+C(1,2,2,1)+C(2,1,1,2))
  print '(3F12.6)', C3
END SUBROUTINE CalcElast



SUBROUTINE CalcPotent(Potential)

! *** Calculation of total potential energy. We assume that the*******
! *** total potential energy can be written as a a superposition of **
! *** central-symmetric forces between two particles.  *************** 

  DOUBLE PRECISION :: R2, R4, R6, R12, Potential, RMin2   

  INTEGER :: PairCnt

  Potential = 0.D0
  DO PairCnt = 1,PairNum
    R2 = PairDist(PairCnt,4)
    IF (R2.LT.RInner*RInner) THEN
      RMin2 = 1/R2
      R4 = RMin2*RMin2
      R6 = R4*RMin2
      R12= R6*R6
      
      Potential = Potential + 4*(R12-R6)
    END IF
  END DO
END SUBROUTINE CalcPotent


SUBROUTINE CalcPairList

! *** Calculation of total potential energy. We assume that the*******
! *** total potential energy can be written as a a superposition of **
! *** central-symmetric forces between two particles.  *************** 

  DOUBLE PRECISION :: R2, D(DIM)

  INTEGER :: I, J, K, PairCnt, CorrCount

  DO PairCnt = 1,PairNum
    I = PairList(PairCnt,1)
    J = PairList(PairCnt,2)
    D = Pos(I,:) -  Pos(J,:)
!    DO K=1, DIM
!      IF (D(K) .GT. VolSize(K)/2) THEN
!        D(K) = D(K) - VolSize(K)
!      ELSE IF (D(K) .LT. -VolSize(K)/2) THEN
!        D(K) = D(K)+VolSize(K)
!      END IF
!    END DO

    R2 = DOT_PRODUCT(D,D)
    PairDist(PairCnt,1:DIM) = D
    PairDist(PairCnt,4) = R2
  END DO
END SUBROUTINE CalcPairList


SUBROUTINE UpdatePairList

! *** Calculation of total potential energy. We assume that the*******
! *** total potential energy can be written as a a superposition of **
! *** central-symmetric forces between two particles.  *************** 

  DOUBLE PRECISION :: R2, D(DIM)

  INTEGER :: I, J, K

  PairNum = 0
  DO I = 1,PartNum
    DO J=I+1,PartNum
      D = Pos(I,:) -  Pos(J,:)
      R2 = DOT_PRODUCT(D,D)
      IF (R2.LT.ROuter*ROuter) THEN
        PairNum = PairNum + 1
        PairList(PairNum,1) = I
        PairList(PairNum,2) = J
      ENDIF
    END DO
  END DO
  IF (PairNum.GT.MaxLen) THEN
    print *, 'not enough memory allocated for pairlist'
  END IF
END SUBROUTINE UpdatePairList


SUBROUTINE CalcTemp(TotEner)

! *** Calculation of total kinetic energy, i.e. the temperature ***
! *** Result is used for rescaling of velocities. *****************


  DOUBLE PRECISION :: TotEner

  INTEGER :: I

  TotEner = SUM(Momentum*Momentum)
  TotEner = TotEner*0.5D0
END SUBROUTINE CalcTemp




SUBROUTINE Rescale (TotEner)

! *** Rescaling velocities to adjust temperature. ****************

  DOUBLE PRECISION :: TotEner, ScaleParam
  INTEGER :: I

  TotEner = TotEner/(PartNum-1)
  ScaleParam = SQRT(Temperature*1.0D0/TotEner)

  Momentum = Momentum*ScaleParam
END SUBROUTINE Rescale 

END MODULE MD2D




