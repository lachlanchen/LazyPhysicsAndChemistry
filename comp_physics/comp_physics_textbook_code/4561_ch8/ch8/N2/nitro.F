      PROGRAM MD
C*****************************************************************
C* This program performs a Molecular Dynamics (MD) simulation    *
C* for nitrogen molecules.                                       *
C* Periodic boundary conditions are assumed. The potential is    *
C* cut off beyond a distance "Rcutoff". Furthermore, the minimum *
C* image convention is used to calculate forces etc.             *
C*                                                               *
C* Global variables are stored in the file "nitro.glob"          *
C*                                                               *
C* Program written by Jos Thijssen	                         *
C* Spring 1999                                                   *
C* This program is described in problem 8.7 of the book          *
C* "Computational Physics" by Jos Thijssen,                      *
C* Cambridge University Press 1999                               *
C*****************************************************************

C**********************Global Variables***************************
       include "nitro.glob"                                       
C*****************************************************************

      LOGICAL Scale
      Call Initialise
      Call InitDraw

C "Scale" controls rescaling of velocities
      
      Scale = .TRUE.
      Call Simulation (Scale, InitTime)
      Scale = .FALSE.
      Call Simulation (Scale, SimTime)

      CALL WriteCorr
      Call WithDraw

      END



      SUBROUTINE InitDraw

C**********************Global Variables***************************
       include "nitro.glob"                                       
C*****************************************************************
      DOUBLE PRECISION WinBor

#ifdef Plot
      CALL InitPlot ("lightblue", 800, 800, "mdplot", 1)
      WinBor = VolSize+0.1D0
      CALL Framing (-0.1D0, -0.1D0, WinBor, WinBor)
      CALL PutStartButton()
      CALL PutStopButton()
#endif
      END

      SUBROUTINE WithDraw

#ifdef Plot
      CALL EndPlot
#endif

      END
 


      SUBROUTINE WriteCorr

C**********************Global Variables***************************
       include "nitro.glob"                                       
C*****************************************************************

      INTEGER I
   
      DO I=1, 250
        WRITE (11, '(I4, F12.6)') I, dble(CorrArray(I))/dble(I*I)
      ENDDO
      CLOSE(2)
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(7)
      END 
     


      SUBROUTINE Initialise

      CALL InitParameters
      CALL InitPositions
      CALL InitMomenta
      CALL OpenFiles

      END


      SUBROUTINE OpenFiles

C**********************Global Variables***************************
       include "nitro.glob"                                       
C*****************************************************************

      OPEN (2, file='potential')
      OPEN (12, file='temperature')
      OPEN (13, file='cmtemp')
      OPEN (7, file='virial')
      OPEN (11, file='correl')
      write (12,*) PartNum
      write (13,*) PartNum
      write (2,*) PartNum
      write (7,*) PartNum
      write (11,*) PartNum
      END

      SUBROUTINE InitParameters

C************  User input ****************************************

C**********************Global Variables***************************
      include  "nitro.glob"                                            
C*****************************************************************

      INTEGER I, Num

      OPEN (UNIT=8, FILE='nitro.in')

      WRITE (6,*) 'Give temperature (in reduced units)'
      READ (8,*) Temperature

      WRITE (6,*) 'Give N-N distance (in red. units)'
      READ (8,*) Sep

      WRITE (6,*) 'Give time step for integration (red. units)'
      READ (8,*) TimeStep

      WRITE (6,*) 'Give equilibration time (red. units)'
      READ (8,*) InitTime

      WRITE (6,*) 'Give nr. of time steps between rescaling'
      READ (8,*) ScaleTime

      WRITE (6,*) 'Give simulation time (red. units)'
      READ (8,*) SimTime

      WRITE (6,*) 'Give density (red. units)'
      READ (8,*) Density

      WRITE (6,*) 'Geef de rng initialisator'
      READ (8,*) num

      WRITE (6,*) 'Geef interval tussen displays'
      READ (8,*) DispInt

C****** Density and particle number define the volume of the system **
      Volume = DBLE(0.5D0*PartNum/Density)
      CorrStep = 0.02
      DO I=1, MaxLen
        CorrArray(I) = 0.D0
      ENDDO

      CALL InitRand(num)

      RInner = 2.5D0
      ROuter = 3.0D0
      

      END


      SUBROUTINE GenDir (DirVec)
     
C**********************Global Variables***************************
      include   "nitro.glob"
C*****************************************************************

      DOUBLE PRECISION Phi, CosTheta, SinTheta, DirVec(3), RealRand, PI

      PI = 4.d0*datan(1.d0)
      Phi = RealRand()*2.d0*PI
      CosTheta = 1.d0-2.d0*RealRand()
      SinTheta = SQRT(1.d0-CosTheta*CosTheta)
      DirVec(1) = SinTheta*COS(Phi)
      DirVec(2) = SinTheta*SIN(Phi)
      DirVec(3) = CosTheta
      END


      SUBROUTINE InitPositions

C *** Positions are stored on a regular fcc lattice. Lattice constant* 
C *** is adjusted to volume size to fill the volume homogeneously ****

C**********************Global Variables***************************
      include   "nitro.glob"
C*****************************************************************

      INTEGER LinCell, IX, IY, IZ, Teller
      DOUBLE PRECISION  LattConst, Third, DirVec(3)


C*******Calculate Volume Size from Volume ************************
C*******LinCell is number of cells along one side ****************
      Third = 1.D0/3.D0
      LinCell = NINT ((DBLE(PartNum)/8)**Third)
      VolSize = Volume**Third
      XVolSize = VolSize
      LattConst = VolSize/LinCell

      Teller = 0
      DO IX = 0, LinCell - 1
        DO IY = 0, LinCell - 1
          DO IZ = 0, LinCell - 1
            Teller = Teller + 1
            CALL GenDir(DirVec)
            OQx(Teller) = (Ix+0.25D0)*LattConst
            OQy(Teller) = (Iy+0.25D0)*LattConst
            OQz(Teller) = (Iz+0.25D0)*LattConst
            Teller = Teller+1
            OQx(Teller) = (Ix+0.25D0)*LattConst+DirVec(1)*Sep
            OQy(Teller) = (Iy+0.25D0)*LattConst+DirVec(2)*Sep
            OQz(Teller) = (Iz+0.25D0)*LattConst+DirVec(3)*Sep
            CALL GenDir(DirVec)
            Teller = Teller + 1
            OQx(Teller) = (Ix+0.75D0)*LattConst
            OQy(Teller) = (Iy+0.75D0)*LattConst
            OQz(Teller) = (Iz+0.25D0)*LattConst
            Teller = Teller+1
            OQx(Teller) = (Ix+0.75D0)*LattConst+DirVec(1)*Sep
            OQy(Teller) = (Iy+0.75D0)*LattConst+DirVec(2)*Sep
            OQz(Teller) = (Iz+0.25D0)*LattConst+DirVec(3)*Sep
            CALL GenDir(DirVec)
            Teller = Teller + 1
            OQx(Teller) = (Ix+0.75D0)*LattConst
            OQy(Teller) = (Iy+0.25D0)*LattConst
            OQz(Teller) = (Iz+0.75D0)*LattConst
            Teller = Teller+1
            OQx(Teller) = (Ix+0.75D0)*LattConst+DirVec(1)*Sep
            OQy(Teller) = (Iy+0.25D0)*LattConst+DirVec(2)*Sep
            OQz(Teller) = (Iz+0.75D0)*LattConst+DirVec(3)*Sep
            CALL GenDir(DirVec)
            Teller = Teller + 1
            OQx(Teller) = (Ix+0.25D0)*LattConst
            OQy(Teller) = (Iy+0.75D0)*LattConst
            OQz(Teller) = (Iz+0.75D0)*LattConst
            Teller = Teller+1
            OQx(Teller) = (Ix+0.25D0)*LattConst+DirVec(1)*Sep
            OQy(Teller) = (Iy+0.75D0)*LattConst+DirVec(2)*Sep
            OQz(Teller) = (Iz+0.75D0)*LattConst+DirVec(3)*Sep
            ENDDO
         ENDDO
      ENDDO
      END



      SUBROUTINE InitMomenta

C     Initialise momenta of particles. All velocity components
C     are drawn from random generator with Gaussian distribution

C**********************Global Variables***************************
      include  "nitro.glob"
C*****************************************************************

      INTEGER I

      DOUBLE PRECISION R1, R2, TotEner, Vx, Vy, Vz

C *****Assign initial velocities to all particles*****************
C *****TotEner is used for rescaling the velocities***************
      TotEner = 0
      DO I = 1, PartNum, 2
         Call ExpRand(R1, R2)
         Px(I) = R1
         Py(I) = R2
         Call ExpRand(R1, R2)
         Pz(I) = R1
         Px(I+1) = R2
         Call ExpRand(R1, R2)
         Py(I+1) = R1
         Pz(I+1) = R2
      ENDDO


C *****Set Centre of Mass velocity equal to zero*******************

      Vx = 0
      Vy = 0
      Vz = 0
      DO I=1,PartNum
        Vx = Vx + Px(I)
        Vy = Vy + Py(I)
        Vz = Vz + Pz(I)
      ENDDO
      Vx = Vx/PartNum
      Vy = Vy/PartNum
      Vz = Vz/PartNum

      DO I=1,PartNum
         Px(I) = Px(I) - Vx
         Py(I) = Py(I) - Vy
         Pz(I) = Pz(I) - Vz
      ENDDO

C ***** Rescale velocities to the right temperature ***************
      Call CalcTemp (TotEner)
      Call Rescale (TotEner)
      DO I=1, PartNum
        Qx(I) = OQx(I) + TimeStep*Px(I)
        Qy(I) = OQy(I) + TimeStep*Py(I)
        Qz(I) = OQz(I) + TimeStep*Pz(I)
      ENDDO
      END






      SUBROUTINE Simulation (Scale, Time)

C *******************************************************************
C ****  Subroutine in which the actual simulation is performed. *****
C **** IF "Scale" is .TRUE., velocities are regularly rescaled ******
C **** "Time" defines simulation time. ******************************
C **** Forces of previous integration step are stored in OldForceX **
C **** etc. ********************************************************* 
C *******************************************************************

C**********************Global Variables***************************
      include   "nitro.glob"
C*****************************************************************

      DOUBLE PRECISION Time,
     .      ForceX(PartNum), ForceY(PartNum), ForceZ(PartNum),
     .      Fac2, TotEner, Virial



      INTEGER StepNum, Step
      LOGICAL Scale, Special, VirCalc

      CALL UpdatePairList
      CALL CalcPairList(Scale)
      CALL CalcForce (ForceX, ForceY, ForceZ, .false.,
     .                Virial)

      StepNum = Time/TimeStep
      Fac2 = 0.5D0*TimeStep*TimeStep
      HalfInvTSt = 0.5D0/TimeStep

      DO Step = 1, StepNum
         IF (DBLE(Step)/ScaleTime .EQ. NINT(DBLE(Step)/ScaleTime)) THEN
            Special = .TRUE.
         ELSE 
            Special = .FALSE.
         ENDIF
         VirCalc =(.NOT.(Scale))
         CALL Integrate(ForceX, ForceY, ForceZ, Fac2, Scale,
     .                  VirCalc, Virial)
            
         IF (DBLE(Step)/DispInt .EQ. NINT(DBLE(Step)/DispInt)) THEN
           CALL DrawConfig (PartNum, VolSize, Qx, Qy, Qz)
         ENDIF
         IF (.not.(Scale)) THEN
           CALL OutputMD(Step, Scale, Virial)
         ENDIF
         IF (Special) THEN
           CALL UpdatePairList
         ENDIF
         IF (Scale .AND. Special) THEN
C *** Rescaling velocities to adjust temperature *******************
            CALL CalcTemp(TotEner)
            CALL Rescale(TotEner)
         ENDIF
      ENDDO

      END



      SUBROUTINE OutputMD(Step, Scale, Virial)

C**********************Global Variables***************************
      include   "nitro.glob"
C*****************************************************************

      INTEGER Step

      LOGICAl Scale

      DOUBLE PRECISION TotEner, Potential, Virial, CMTemp

      CALL CalcTemp (TotEner)
      CALL CalcCMTemp (CMTemp)
      CALL CalcPotent (Potential)
         WRITE (12,'(2F12.5)') Step*TimeStep+InitTime, 
     .            TotEner*0.8D0/(PartNum-1)
         WRITE (13,'(2F12.5)') Step*TimeStep+InitTime, 
     .            CMTemp/(PartNum-1)*1.3333333333d0
         WRITE (2,'(2F12.5)') Step*TimeStep+InitTime,  
     .                              Potential
         WRITE (7, '(2F12.5)') Step*TimeStep+InitTime,  
     .                        Virial
      IF (DBLE(Step)/DispInt .EQ. NINT(DBLE(Step)/DispInt)) THEN
        ENDFILE 7
        ENDFILE 2
        ENDFILE 12
      ENDIF
      END


      SUBROUTINE DrawConfig (PartNum, VolSize, Qx, Qy, Qz)

      INTEGER PartNum, I
      DOUBLE PRECISION VolSize, QX(PartNum), QY(PartNum), QZ(PartNum)
#ifdef Plot
      DO I=1, PartNum/2
        dist = (qx(2*i-1)-qx(2*I))**2 + 
     .         (qy(2*i-1)-qy(2*I))**2 +
     .         (qz(2*i-1)-qz(2*I))**2
	CALL Draw (Qx(2*I-1), Qy(2*I-1), Qx(2*I), Qy(2*I))
      ENDDO
#endif
      END


      DOUBLE PRECISION FUNCTION Norm(X, Y, Z)
      IMPLICIT NONE
      DOUBLE PRECISION X, Y, Z
      Norm = X*X+Y*Y+Z*Z
      END
      
      DOUBLE PRECISION FUNCTION InProd(X1, Y1, Z1, X2, Y2, Z2)
      IMPLICIT NONE
      DOUBLE PRECISION X1, Y1, Z1, X2, Y2, Z2
      InProd = X1*X2+Y1*Y2+Z1*Z2
      END
 

      DOUBLE PRECISION FUNCTION Lambda(A, B, AB)

C**********************Global Variables***************************
      include   "nitro.glob"
C*****************************************************************

      DOUBLE PRECISION A, B, AB, L1, L2, D
      D = SQRT(AB*AB-A*B)

      L1 = (-AB+D)/B*0.5D0
      L2 = (-AB-D)/B*0.5D0
      IF (ABS(L1).LT. ABS(L2)) THEN
        Lambda = L1
      ELSE
        Lambda = L2
      ENDIF
      END



      SUBROUTINE Integrate (ForceX, ForceY, ForceZ,
     .                      Fac2, Scale, VirClc, Virial)


C *** Integration of equations of motion using Verlet algorithm ***

C**********************Global Variables***************************
      include   "nitro.glob"
C*****************************************************************

      DOUBLE PRECISION ForceX(PartNum), ForceY(PartNum), 
     .       ForceZ(PartNum),
     .       Fac2, Virial, A1x, A1y, A1z, A2x, A2y, A2z,
     .       A12x, A12y, A12z, SqTSt, AvQx, AvQy, AvQz, 
     .       B12x, B12y, B12z, A, B, AB, InProd, Norm, Lambda, Param

      INTEGER I
      LOGICAL Scale, VirClc

      SqTSt = TimeStep**2
      DO I=1,PartNum/2
         A1x = 2*Qx(2*I-1)-OQx(2*I-1)+SqTSt*ForceX(2*I-1)
         A2x = 2*Qx(2*I)  -OQx(2*I)  +SqTSt*ForceX(2*I)
         A1y = 2*Qy(2*I-1)-OQy(2*I-1)+SqTSt*ForceY(2*I-1)
         A2y = 2*Qy(2*I)  -OQy(2*I)  +SqTSt*ForceY(2*I)
         A1z = 2*Qz(2*I-1)-OQz(2*I-1)+SqTSt*ForceZ(2*I-1)
         A2z = 2*Qz(2*I)  -OQz(2*I)  +SqTSt*ForceZ(2*I)
         A12x = A1x-A2x
         A12y = A1y-A2y
         A12z = A1z-A2z
         B12x = (Qx(2*I-1)-Qx(2*I))*SqTSt
         B12y = (Qy(2*I-1)-Qy(2*I))*SqTSt
         B12z = (Qz(2*I-1)-Qz(2*I))*SqTSt
         A = Norm(A12X, A12Y, A12Z)-Sep*Sep
         B = Norm(B12X, B12Y, B12Z)
         AB = InProd(A12X, A12Y, A12Z, B12X, B12Y, B12Z)
         Param = Lambda(A, B, AB)
C In the next few lines, the "old" values of the positions are temporarily
C stored in the P's in order to be able to calculate the latter below.
         Px(2*I-1)  = OQx(2*I-1)
         Py(2*I-1)  = OQy(2*I-1)
         Pz(2*I-1)  = OQz(2*I-1)
         Px(2*I)    = OQx(2*I)
         Py(2*I)    = OQy(2*I)
         Pz(2*I)    = OQz(2*I)
         OQx(2*I-1) = Qx(2*I-1)
         OQy(2*I-1) = Qy(2*I-1)
         OQz(2*I-1) = Qz(2*I-1)
         OQx(2*I)   = Qx(2*I)
         OQy(2*I)   = Qy(2*I)
         OQz(2*I)   = Qz(2*I)
         Qx(2*I-1)  = A1x+Param*B12x
         Qy(2*I-1)  = A1y+Param*B12y
         Qz(2*I-1)  = A1z+Param*B12z
         Qx(2*I)    = A2x-Param*B12x
         Qy(2*I)    = A2y-Param*B12y
         Qz(2*I)    = A2z-Param*B12z
         Px(2*I-1)  = (Qx(2*I-1)-Px(2*I-1))*HalfInvTSt
         Py(2*I-1)  = (Qy(2*I-1)-Py(2*I-1))*HalfInvTSt
         Pz(2*I-1)  = (Qz(2*I-1)-Pz(2*I-1))*HalfInvTSt
         Px(2*I)    = (Qx(2*I)-Px(2*I))*HalfInvTSt
         Py(2*I)    = (Qy(2*I)-Py(2*I))*HalfInvTSt
         Pz(2*I)    = (Qz(2*I)-Pz(2*I))*HalfInvTSt
         AvQx = 0.5D0*(Qx(2*I-1)+Qx(2*I))
         AvQy = 0.5D0*(Qy(2*I-1)+Qy(2*I))
         AvQz = 0.5D0*(Qz(2*I-1)+Qz(2*I))
         IF (AvQx .LT. 0.D0) THEN
            Qx(2*I) = Qx(2*I) + VolSize
            Qx(2*I-1) = Qx(2*I-1) + VolSize
            OQx(2*I) = OQx(2*I) + VolSize
            OQx(2*I-1) = OQx(2*I-1) + VolSize
         ELSE
            IF (AvQx .GT. VolSize) THEN
              Qx(2*I) = Qx(2*I) - VolSize
              Qx(2*I-1) = Qx(2*I-1) - VolSize
              OQx(2*I) = OQx(2*I) - VolSize
              OQx(2*I-1) = OQx(2*I-1) - VolSize
            ENDIF
         ENDIF
         IF (AvQy .LT. 0.D0) THEN
            Qy(2*I) = Qy(2*I) + VolSize
            Qy(2*I-1) = Qy(2*I-1) + VolSize
            OQy(2*I) = OQy(2*I) + VolSize
            OQy(2*I-1) = OQy(2*I-1) + VolSize
         ELSE
            IF (AvQy .GT. VolSize) THEN
              Qy(2*I) = Qy(2*I) - VolSize
              Qy(2*I-1) = Qy(2*I-1) - VolSize
              OQy(2*I) = OQy(2*I) - VolSize
              OQy(2*I-1) = OQy(2*I-1) - VolSize
            ENDIF
         ENDIF
         IF (AvQz .LT. 0.D0) THEN
            Qz(2*I) = Qz(2*I) + VolSize
            Qz(2*I-1) = Qz(2*I-1) + VolSize
            OQz(2*I) = OQz(2*I) + VolSize
            OQz(2*I-1) = OQz(2*I-1) + VolSize
         ELSE
            IF (AvQz .GT. VolSize) THEN
              Qz(2*I) = Qz(2*I) - VolSize
              Qz(2*I-1) = Qz(2*I-1) - VolSize
              OQz(2*I) = OQz(2*I) - VolSize
              OQz(2*I-1) = OQz(2*I-1) - VolSize
            ENDIF
         ENDIF
      ENDDO
      CALL CalcPairList(Scale)
      CALL CalcForce (ForceX,ForceY,ForceZ, VirClc, Virial)
      END





      SUBROUTINE CalcForce(ForceX, ForceY, ForceZ, VirClc, Virial)

C *** Calculation of forces. We assume that the forces are a super- **
C *** position of central-symmetric forces between two particles.   ** 

C**********************Global Variables***************************
      include   "nitro.glob"
C*****************************************************************

      DOUBLE PRECISION ForceX(PartNum), ForceY(PartNum), 
     .       ForceZ(PartNum),
     .       R2, R4, R8, R14, ForceConst, Sx, Sy, Sz, 
     .       RMin2, Virial

      LOGICAL VirClc

      INTEGER I, J, PairCnt

      DO I=1, PartNum
        ForceX(I) = -.0D0
        ForceY(I) = -.0D0
        ForceZ(I) = -.0D0
      ENDDO
      Virial = 0.D0
      DO PairCnt = 1,PairNum
        I = PairList(PairCnt,1)
        J = PairList(PairCnt,2)
        R2 = PairDist(PairCnt,4)
        IF (R2.LT.Rinner*RInner) THEN
          Sx = PairDist(PairCnt,1)
          Sy = PairDist(PairCnt,2)
          Sz = PairDist(PairCnt,3)
          RMin2 = 1/R2
          R4 = RMin2*RMin2
          R8 = R4*R4
          R14= R8*R4*RMin2
          ForceConst = 12.d0*(R14-R8)
          ForceX(I) = ForceX(I) + Sx*ForceConst
          ForceY(I) = ForceY(I) + Sy*ForceConst
          ForceZ(I) = ForceZ(I) + Sz*ForceConst

          ForceX(J) = ForceX(J) - Sx*ForceConst
          ForceY(J) = ForceY(J) - Sy*ForceConst
          ForceZ(J) = ForceZ(J) - Sz*ForceConst
          Virial = Virial + ForceConst*R2
        ENDIF
      ENDDO
      
      END


      LOGICAL FUNCTION SameMolec(I, J)

C**********************Global Variables***************************
      include   "nitro.glob" 
C*****************************************************************

      INTEGER I, J

      SameMolec = ((MOD(J,2).EQ.0) .AND. (J.EQ.I+1))

      END



      SUBROUTINE CalcPotent(Potential)

C *** Calculation of total potential energy. We assume that the*******
C *** total potential energy can be written as a a superposition of **
C *** central-symmetric forces between two particles.  *************** 

C**********************Global Variables***************************
      include   "nitro.glob" 
C*****************************************************************

      DOUBLE PRECISION R2, R4, R6, R12, Potential,
     .       RMin2   

      INTEGER PairCnt

      Potential = 0.D0
      DO PairCnt = 1,PairNum
        R2 = PairDist(PairCnt,4)
        IF (R2.LT.RInner*RInner) THEN
          RMin2 = 1/R2
          R4 = RMin2*RMin2
          R6 = R4*RMin2
          R12= R6*R6
         
          Potential = Potential + (R12-2.d0*R6)
        ENDIF
      ENDDO
      END

      SUBROUTINE CalcPairList(Scale)

C *** Calculation of total potential energy. We assume that the*******
C *** total potential energy can be written as a a superposition of **
C *** central-symmetric forces between two particles.  *************** 

C**********************Global Variables***************************
      include   "nitro.glob" 
C*****************************************************************

      DOUBLE PRECISION R2, Dx, Dy, Dz

      INTEGER I, J, PairCnt, CorrCount

      LOGICAL Scale

      DO PairCnt = 1,PairNum
        I = PairList(PairCnt,1)
        J = PairList(PairCnt,2)
        Dx = Qx(I) -  Qx(J)
        IF (Dx .GT. XVolSize/2) THEN
            Dx = Dx - XVolSize
        ENDIF
        IF (Dx .LT. -XVolSize/2) THEN
            Dx = Dx+XVolSize
        ENDIF

        Dy = Qy(I) -  Qy(J)
        IF (Dy .GT. VolSize/2) THEN
            Dy = Dy - VolSize
        ENDIF
        IF (Dy .LT. -VolSize/2) THEN
            Dy = Dy+VolSize
        ENDIF
        Dz = Qz(I) -  Qz(J)
        IF (Dz .GT. VolSize/2) THEN
          Dz = Dz - VolSize
        ENDIF
        IF (Dz .LT. -VolSize/2) THEN
          Dz = Dz+VolSize
        ENDIF

        R2 = Dx*Dx+Dy*Dy+Dz*Dz
        PairDist(PairCnt,1) = Dx
        PairDist(PairCnt,2) = Dy
        PairDist(PairCnt,3) = Dz
        PairDist(PairCnt,4) = R2
        IF (.not.(Scale)) THEN
          CorrCount = nint(SQRT(R2)/CorrStep)+1
          CorrArray(CorrCount) = CorrArray(CorrCount)+1
        ENDIF
      ENDDO
      END


      SUBROUTINE UpdatePairList

C *** Calculation of total potential energy. We assume that the*******
C *** total potential energy can be written as a a superposition of **
C *** central-symmetric forces between two particles.  *************** 

C**********************Global Variables***************************
      include   "nitro.glob" 
C*****************************************************************

      DOUBLE PRECISION R2, Dx, Dy, Dz

      INTEGER I, J, StartJ, misnum

      PairNum = 0
      MisNum = 0
      DO I = 1,PartNum
        IF (MOD(I,2) .EQ.1) THEN
          StartJ = I+2
        ELSE
          StartJ = I+1
        ENDIF
        DO J=StartJ,PartNum
           Dx = Qx(I) -  Qx(J)
           IF (Dx .GT. XVolSize/2) THEN
               Dx = Dx - XVolSize
           ENDIF
           IF (Dx .LT. -XVolSize/2) THEN
               Dx = Dx+XVolSize
           ENDIF

           Dy = Qy(I) -  Qy(J)
           IF (Dy .GT. VolSize/2) THEN
               Dy = Dy - VolSize
           ENDIF
           IF (Dy .LT. -VolSize/2) THEN
               Dy = Dy+VolSize
           ENDIF
           Dz = Qz(I) -  Qz(J)
           IF (Dz .GT. VolSize/2) THEN
               Dz = Dz - VolSize
           ENDIF
           IF (Dz .LT. -VolSize/2) THEN
               Dz = Dz+VolSize
           ENDIF

           R2 = Dx*Dx+Dy*Dy+Dz*Dz
           IF (R2.LT.ROuter*ROuter) THEN
             PairNum = PairNum + 1
             PairList(PairNum,1) = I
             PairList(PairNum,2) = J
           else
             misnum = misnum+1
           ENDIF
        ENDDO
      ENDDO
      IF (PairNum.GT.MaxLen) THEN
        print *, 'not enough memory allocated for pairlist'
      ENDIF
      END


      SUBROUTINE CalcTemp(TotEner)

C *** Calculation of total kinetic energy, i.e. the temperature ***
C *** Result is used for rescaling of velocities. *****************

C**********************Global Variables***************************
      include   "nitro.glob"
C*****************************************************************

      DOUBLE PRECISION TotEner

      INTEGER I

      TotEner = 0.D0
      DO I=1, PartNum
         TotEner = TotEner + Px(I)*Px(I) + Py(I)*Py(I)+
     .                       Pz(I)*Pz(I)
      ENDDO
      TotEner = TotEner*0.5D0
      END

      SUBROUTINE CalcCMTemp(CMEner)

C *** Calculation of total kinetic energy, i.e. the temperature ***
C *** Result is used for rescaling of velocities. *****************

C**********************Global Variables***************************
      include   "nitro.glob"
C*****************************************************************

      DOUBLE PRECISION CMEner, PCMx, PCMy, PCMz

      INTEGER I

      CMEner = 0.D0
      DO I=1, PartNum/2
         PCMx = 0.5D0*(Px(2*I-1)+Px(2*I))
         PCMy = 0.5D0*(Py(2*I-1)+Py(2*I))
         PCMz = 0.5D0*(Pz(2*I-1)+Pz(2*I))
         CMEner = CMEner + PCMx*PCMx+PCMy*PCMy+PCMz*PCMz
      ENDDO
C No factor 0.5 because the mass of a molecule is twice that of an atom

      END




      SUBROUTINE Rescale (TotEner)

C *** Rescaling velocities to adjust temperature. ****************

C**********************Global Variables***************************
C*                                                               *
      include  "nitro.glob"
C*                                                               *
C*                 Stored in file: "nitro.glob"                   *
C*  Px, Py, Pz, Qx, Qy, Qz : DOUBLE PRECISION arrays of size partnum       *
C*  PartNum, Temperature, Volume, InitTime, ScaleTime, TimeStep, *
C*  SimTime                                                      *
C*****************************************************************

      DOUBLE PRECISION TotEner, ScaleParam
      INTEGER I

C      return   

      TotEner = TotEner/(PartNum-1)
      ScaleParam = SQRT(Temperature*1.25D0/TotEner)

      DO I=1,PartNum
         Px(I) = Px(I)*ScaleParam
         Py(I) = Py(I)*ScaleParam
         Pz(I) = Pz(I)*ScaleParam
      ENDDO
      END






