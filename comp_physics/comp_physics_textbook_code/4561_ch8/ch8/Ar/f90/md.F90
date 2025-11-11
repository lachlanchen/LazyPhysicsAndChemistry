PROGRAM MD
!*****************************************************************
!* This program performs a Molecular Dynamics (MD) simulation    *
!* for a monatomic material.                                     *
!* Periodic boundary conditions are assumed. The potential is    *
!* cut off beyond a distance "Rcutoff". Furthermore, the minimum *
!* image convention is used to calculate forces etc.             *
!*                                                               *
!*                                                               *
!* Program written by Jos Thijssen	                         *
!* Program written 2004-2007                                     *
!* For input parameters, see routine 'InitParameters' below      *
!*****************************************************************
  IMPLICIT NONE

  INTEGER, PARAMETER :: PartNum = 864, MaxLen=PartNum*(PartNum-1)/2, Dim=3, MaxCorr=250

  DOUBLE PRECISION :: Pos(PartNum, Dim), S, &
            Momentum(PartNum, Dim), &
            Temperature, Volume, InitTime, SimTime, TimeStep, &
            VolSize, Density, MaxRange, CorrStep, &
            PairDist(MaxLen,4), ROuter, RInner, Sep 

  INTEGER :: ScaleTime, CorrArray(MaxCorr), PairList(MaxLen, 2), &
           PairNum, DispInt

  LOGICAL :: Scale

  Call Initialise

  Scale = .TRUE.
  Call Simulation (Scale, InitTime)
  Scale = .FALSE.
  Call Simulation (Scale, SimTime) 

  CALL FinalWrite

CONTAINS


SUBROUTINE FinalWrite

 INTEGER :: I

! Write correlation function
 DO I=1, MaxCorr
   WRITE (11, '(I4, F12.6)') I, dble(CorrArray(I))/dble(I*I)
 ENDDO
! Close output files
 CLOSE(2)
 CLOSE(11)
 CLOSE(12)
 CLOSE(7)
#ifdef Plot
 CALL EndPlot()
#endif
END SUBROUTINE FinalWrite


SUBROUTINE Initialise
  CALL InitParameters
  CALL InitPositions
  CALL InitMomenta
  CALL OpenFiles
#ifdef Plot
  CALL InitPlot('lightblue', 700,700, 'out.ps', 1)
  CALL Framing(-0.05D0*VolSize, -0.05*VolSize, 1.05D0*VolSize, 1.05*VolSize)
  CALL PutStopButton
#endif
END SUBROUTINE Initialise


SUBROUTINE OpenFiles
 OPEN (2, file='potential')
 OPEN (12, file='totener')
 OPEN (8, file='temperature')
 OPEN (7, file='virial')
 OPEN (11, file='correl')
END SUBROUTINE OpenFiles



SUBROUTINE InitParameters

  INTEGER :: I, Num

  OPEN (UNIT=8, FILE='md.in')

  READ (8,*) Temperature
  READ (8,*) TimeStep
  READ (8,*) InitTime
  READ (8,*) ScaleTime
  READ (8,*) SimTime
  READ (8,*) Density
  READ (8,*) num
  READ (8,*) DispInt

!****** Density and particle number define the volume of the system **
  Volume = DBLE(PartNum/Density)
  WRITE (6,*) 'Volume = ', Volume
  CorrStep = 0.02D0
  DO I=1, MaxCorr
    CorrArray(I) = 0
  END DO
  CALL InitRand(num)
  RInner = 3.0D0
  ROuter = 3.7D0
END SUBROUTINE InitParameters



SUBROUTINE InitPositions

! *** Positions are stored on a regular fcc lattice. Lattice constant* 
! *** is adjusted to volume size to fill the volume homogeneously ****

INTEGER :: LinCell, IX, IY, IZ, Counter
DOUBLE PRECISION :: LattConst, Third

!*******Calculate Volume Size from Volume ************************
!*******LinCell is number of cells along one side ****************
  Third = 1.D0/3.D0
  WRITE (6,*) Volume, PartNum
  LinCell = NINT ((DBLE(PartNum)/4)**Third)
  WRITE (6,*) 'LinCell = ', LinCell
  VolSize = Volume**Third
  LattConst = VolSize/LinCell
  WRITE (6,*) 'LattConst = ', LattConst, LinCell
  Counter = 0
  DO IX = 0, LinCell - 1
    DO IY = 0, LinCell - 1
      DO IZ = 0, LinCell - 1
        Counter = Counter + 1
        Pos(Counter,:) = (/(Ix+0.25D0)*LattConst, (Iy+0.25D0)*LattConst, &
                           (Iz+0.25D0)*LattConst/)
        Counter = Counter + 1
        Pos(Counter,:) = (/(Ix+0.75D0)*LattConst, (Iy+0.75D0)*LattConst, &
                           (Iz+0.25D0)*LattConst/)
        Counter = Counter + 1
        Pos(Counter,:) = (/(Ix+0.75D0)*LattConst, (Iy+0.25D0)*LattConst, &
                           (Iz+0.75D0)*LattConst/)
        Counter = Counter + 1
        Pos(Counter,:) = (/(Ix+0.25D0)*LattConst, (Iy+0.75D0)*LattConst, &
                           (Iz+0.75D0)*LattConst/)
      END DO
    END DO
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
  DO I = 1, PartNum, 2
    Call ExpRand(R1, R2)
    Momentum(I,1:2) = (/R1, R2/)
    Call ExpRand(R1, R2)
    Momentum(I,3) = R1
    Momentum(I+1,1) = R2
    Call ExpRand(R1, R2)
    Momentum(I+1,2:3) = (/R1, R2/)
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

#ifdef Plot
SUBROUTINE PlotConf
! Plot the configuration graphically, projected onto XY plane
  INTEGER :: I
  DO I=1, PartNum
    CALL SetPoint(Pos(I,1), Pos(I,2))
  END DO
END SUBROUTINE PlotConf
#endif

SUBROUTINE Simulation (Scale, Time)

! *******************************************************************
! ****  Subroutine in which the actual simulation is performed. *****
! **** IF "Scale" is .TRUE., velocities are regularly rescaled ******
! **** "Time" defines simulation time. ******************************
! **** Forces of previous integration step are stored in OldForceX **
! **** etc. ********************************************************* 
! *******************************************************************

  DOUBLE PRECISION :: Time, Force(PartNum,DIM), &
                      KinEner, Virial

  INTEGER :: StepNum, Step
  LOGICAL :: Scale

  CALL UpdatePairList
  CALL CalcPairList(Scale)
  CALL CalcForce (Force, Virial)

  StepNum = INT(Time/TimeStep)
  WRITE (6,*) StepNum
  DO Step = 1, StepNum
! Special is a Boolean variable which determines whether rescaling
! and/or calculation of virial takes place during this particular step. 
    CALL Integrate(Force, Scale, Virial)
    IF(MODULO(Step,ScaleTime)==0.AND.Scale) THEN
      CALL CalcTemp(KinEner)
      CALL Rescale(KinEner)
    END IF
    IF (MODULO(Step,DispInt)==0) THEN
      WRITE (6,*) Step
      CALL UpdatePairList
#ifdef Plot
      CALL PlotConf()
#endif
      IF (.NOT.Scale) THEN
        CALL OutputMD(Step, Virial)
      END IF
    END IF
  END DO
END SUBROUTINE Simulation



SUBROUTINE OutputMD(Step, Virial)

  INTEGER :: Step

  DOUBLE PRECISION :: Potential, KinEner, Virial

  CALL CalcTemp (KinEner)
  CALL CalcPotent (Potential)
  WRITE (12,'(I8,F12.5)') Step, KinEner+Potential
  WRITE (8,'(F12.5)') KinEner
  WRITE (2, '(3F12.5)') Potential
  WRITE (7, '(F12.5)') Virial
#ifdef Plot
!  CALL SetNamedBackground('lightblue')
#endif
END SUBROUTINE OutputMD


SUBROUTINE Integrate (Force, Scale, Virial)


! *** Integration of equations of motion using velocity-Verlet algorithm ***

  DOUBLE PRECISION :: Force(PartNum, DIM), &
                    Fac2, TotKin, Virial

  INTEGER :: I, K
  LOGICAL :: Scale

  Momentum = Momentum + 0.5D0*TimeStep*Force
  Pos = Pos + TimeStep*Momentum
! Positions are moved back into unit cell if necessary.
  Pos = MODULO(Pos, VolSize)

  CALL CalcPairList(Scale)
  CALL CalcForce (Force, Virial)
  Momentum = Momentum + 0.5D0*TimeStep*Force
  CALL CalcTemp(TotKin)
END SUBROUTINE Integrate 





SUBROUTINE CalcForce(Force, Virial)

! *** Calculation of forces. We assume that the forces are a super- **
! *** position of central-symmetric forces between two particles.   ** 

  DOUBLE PRECISION :: Force(PartNum,Dim),&
                    R2, R4, R8, R14, ForceConst, S(DIM), &
                    RMin2, Virial 

  INTEGER :: I, J, PairCnt

  Force = 0.D0
  Virial = 0.D0
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
      Force(I,:) = Force(I,:) + S*ForceConst

      Force(J,:) = Force(J,:) - S*ForceConst
      Virial = Virial + ForceConst*R2
    END IF
  END DO
END SUBROUTINE CalcForce






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
      Potential = Potential +4*(R12-R6)
    END IF
  END DO
END SUBROUTINE CalcPotent


SUBROUTINE CalcPairList(Scale)

! *** Calculation of total potential energy. We assume that the*******
! *** total potential energy can be written as a a superposition of **
! *** central-symmetric forces between two particles.  *************** 

  DOUBLE PRECISION :: R2, D(DIM)

  INTEGER :: I, J, K, PairCnt, CorrCount

  LOGICAL :: Scale

  DO PairCnt = 1,PairNum
    I = PairList(PairCnt,1)
    J = PairList(PairCnt,2)
    D = Pos(I,:) -  Pos(J,:)
    DO K=1, DIM
      IF (D(K) .GT. VolSize/2) THEN
        D(K) = D(K) - VolSize
      ELSE IF (D(K) .LT. -VolSize/2) THEN
        D(K) = D(K)+VolSize
      END IF
    END DO

    R2 = DOT_PRODUCT(D,D)
    PairDist(PairCnt,1:DIM) = D
    PairDist(PairCnt,4) = R2
    IF (.NOT.Scale) THEN
      CorrCount = nint(SQRT(R2)/CorrStep)+1
      IF (CorrCount<MaxCorr) THEN
        CorrArray(CorrCount) = CorrArray(CorrCount)+1
      END IF
    END IF
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
! The following is not the most compact way of accounting for periodic 
! boundary conditions, but it is probably the most efficient.
      D = Pos(I,:) -  Pos(J,:)
      DO K=1, DIM
        IF (D(K) .GT. VolSize/2) THEN
          D(K) = D(K) - VolSize
        ELSE IF (D(K) .LT. -VolSize/2) THEN
          D(K) = D(K)+VolSize
        END IF
      END DO

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


SUBROUTINE CalcTemp(KinEner)

! *** Calculation of total kinetic energy, i.e. the temperature ***
! *** Result is used for rescaling of velocities. *****************


  DOUBLE PRECISION :: KinEner

  INTEGER :: I

  KinEner = SUM(Momentum*Momentum)
  KinEner = KinEner*0.5D0
END SUBROUTINE CalcTemp




SUBROUTINE Rescale (KinEner)

! *** Rescaling velocities to adjust temperature. ****************

  DOUBLE PRECISION :: KinEner, ScaleParam
  INTEGER :: I

  KinEner = KinEner/(PartNum-1)
  ScaleParam = SQRT(Temperature*1.5D0/KinEner)

  Momentum = Momentum*ScaleParam
END SUBROUTINE Rescale 

END PROGRAM MD




