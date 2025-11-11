      PROGRAM MD
C*****************************************************************
C* This program performs a Molecular Dynamics (MD) simulation    *
C* for a monatomic material.                                     *
C* Periodic boundary conditions are assumed. The potential is    *
C* cut off beyond a distance "Rcutoff". Furthermore, the minimum *
C* image convention is used to calculate forces etc.             *
C*                                                               *
C* Global variables are stored in the file "md.glb"              *
C* Program described in "Computational Physics", J. M. Thijssen  *
C* Section 8.3                                                   *
C*                                                               *
C* Program written by Jos Thijssen	                         *
C* Program written 1991-1999                                     *
C*****************************************************************

C**********************Global Variables***************************
       include "md.glb"                                       
C*****************************************************************

      LOGICAL Scale

      Call Initialise

C "Scale" controls rescaling of velocities
      Scale = .TRUE.
C Yes, rescale velocities
      Call Simulation (Scale, InitTime)
      Scale = .FALSE.
C No, do not rescale
      Call Simulation (Scale, SimTime)

      CALL FinalWrite

      END


 


      SUBROUTINE FinalWrite

C**********************Global Variables***************************
       include "md.glb"                                       
C*****************************************************************

      INTEGER I
   
C Write correlation function
      DO I=1, 250
        WRITE (11, '(I4, F12.6)') I, dble(CorrArray(I))/dble(I*I)
      ENDDO
C Close output files
      CLOSE(2)
      CLOSE(11)
      CLOSE(12)
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
       include "md.glb"                                       
C*****************************************************************

      OPEN (2, file='potential')
      OPEN (12, file='temperature')
      OPEN (7, file='virial')
      OPEN (11, file='correl')
      END



      SUBROUTINE InitParameters

C************  User input ****************************************

C**********************Global Variables***************************
      include  "md.glb"                                            
C*****************************************************************

      INTEGER I, Num

      WRITE (6,1000) PartNum
      OPEN (UNIT=8, FILE='md.in')

      READ (8,*) Temperature
      READ (8,*) TimeStep
      READ (8,*) InitTime
      READ (8,*) ScaleTime
      READ (8,*) SimTime
      READ (8,*) Dichtheid
      READ (8,*) num
      READ (8,*) DispInt

C****** Density and particle number define the volume of the system **
      Volume = DBLE(PartNum/Dichtheid)
      WRITE (6,*) 'Volume = ', Volume
1000  FORMAT ('This is a simulation with ', I3,' particles.')
      CorrStep = 0.02D0
      DO I=1, MaxLen
        CorrArray(I) = 0
      ENDDO

      CALL InitRand(num)

C Radii for Verlet's neighbour list
      RInner = 3.0D0
      ROuter = 3.7D0
      

      END




      SUBROUTINE InitPositions

C *** Positions are stored on a regular fcc lattice. Lattice constant* 
C *** is adjusted to volume size to fill the volume homogeneously ****

C**********************Global Variables***************************
      include   "md.glb"
C*****************************************************************

      INTEGER LinCell, IX, IY, IZ, Counter
      DOUBLE PRECISION  LattConst, Third


C*******Calculate Volume Size from Volume ************************
C*******LinCell is number of cells along one side ****************
      Third = 1.D0/3.D0
      WRITE (6,*) Volume, PartNum
      LinCell = NINT ((DBLE(PartNum)/4)**Third)
      WRITE (6,*) 'LinCell = ', LinCell
      VolSize = Volume**Third
      XVolSize = VolSize
      LattConst = VolSize/LinCell
      WRITE (6,*) 'LattConst = ', LattConst, LinCell

      Counter = 0
      DO IX = 0, LinCell - 1
         DO IY = 0, LinCell - 1
            DO IZ = 0, LinCell - 1
               Counter = Counter + 1
               Qx(Counter) = (Ix+0.25D0)*LattConst
               Qy(Counter) = (Iy+0.25D0)*LattConst
               Qz(Counter) = (Iz+0.25D0)*LattConst
               Counter = Counter + 1
               Qx(Counter) = (Ix+0.75D0)*LattConst
               Qy(Counter) = (Iy+0.75D0)*LattConst
               Qz(Counter) = (Iz+0.25D0)*LattConst
               Counter = Counter + 1
               Qx(Counter) = (Ix+0.75D0)*LattConst
               Qy(Counter) = (Iy+0.25D0)*LattConst
               Qz(Counter) = (Iz+0.75D0)*LattConst
               Counter = Counter + 1
               Qx(Counter) = (Ix+0.25D0)*LattConst
               Qy(Counter) = (Iy+0.75D0)*LattConst
               Qz(Counter) = (Iz+0.75D0)*LattConst
            ENDDO
         ENDDO
      ENDDO
      END



      SUBROUTINE InitMomenta

C     Initialise momenta of particles. All velocity components
C     are drawn from random generator with Gaussian distribution

C**********************Global Variables***************************
      include  "md.glb"
C*****************************************************************

      INTEGER I

      DOUBLE PRECISION R1, R2, TotEner, Vx, Vy, Vz

C *****Assign initial velocities to all particles*****************
C *****TotEner is used for rescaling the velocities***************
      TotEner = 0
      DO I = 1, PartNum, 2
         Call ExpRand(R1, R2)
c         print *, r1, r2
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
      include   "md.glb"
C*****************************************************************

      DOUBLE PRECISION Time,
     .      ForceX(PartNum), ForceY(PartNum), ForceZ(PartNum),
     .      KinEner, Virial

      INTEGER StepNum, Step
      LOGICAL Scale

      CALL UpdatePairList
      CALL CalcPairList(Scale)
      CALL CalcForce (ForceX, ForceY, ForceZ,
     .                Virial)

      StepNum = INT(Time/TimeStep)
      WRITE (6,*) StepNum
      DO Step = 1, StepNum
        CALL Integrate(ForceX, ForceY, ForceZ, Scale, Virial)
        IF ((MOD(Step, ScaleTime).EQ.0).AND.Scale) THEN
          CALL CalcTemp(KinEner)
          CALL Rescale(KinEner)
        END IF
        IF (MOD(Step, DispInt) .EQ.0) THEN
          WRITE (6,*) Step
          CALL UpdatePairList
          IF (.NOT.(Scale)) THEN
             CALL OutputMD(Step, Virial)
          END IF
        ENDIF
      ENDDO

      END



      SUBROUTINE OutputMD(Step, Virial)

C**********************Global Variables***************************
      include   "md.glb"
C*****************************************************************

      INTEGER Step

      DOUBLE PRECISION KinEner, Potential, Virial

      CALL CalcTemp (KinEner)
      CALL CalcPotent (Potential)
      WRITE (12,'(F12.5)') KinEner
      WRITE (2, '(F12.5)') Potential
      WRITE (7, '(F12.5)') Virial
C The following statements ensure that the ouput buffers are regularly written
C to the output files, so that stopping the program does not cause loss of all
C data
      IF (DBLE(Step)/DispInt .EQ. NINT(DBLE(Step)/DispInt)) THEN
        ENDFILE 7
        ENDFILE 2
        ENDFILE 12
      END IF
      END



      SUBROUTINE Integrate (ForceX, ForceY, ForceZ,
     .                      Scale, Virial)


C *** Integration of equations of motion using velocity-Verlet algorithm ***

C**********************Global Variables***************************
      include   "md.glb"
C*****************************************************************

      DOUBLE PRECISION ForceX(PartNum), ForceY(PartNum), 
     .       ForceZ(PartNum),
     .       Fac2, TotKin, Virial

      INTEGER I
      LOGICAL Scale

      DO I=1,PartNum   
         Px(I) = Px(I) + 0.5D0*TimeStep*ForceX(I)
         Py(I) = Py(I) + 0.5D0*TimeStep*ForceY(I)
         Pz(I) = Pz(I) + 0.5D0*TimeStep*ForceZ(I)
      ENDDO
      DO I=1,PartNum
         Qx(I) = Qx(I) + TimeStep*Px(I)
         Qy(I) = Qy(I) + TimeStep*Py(I)
         Qz(I) = Qz(I) + TimeStep*Pz(I)
C Positions are moved back into unit cell if necessary.
C Note that this can be done using the code Qx(I) = Qx(I) - INT(Qx/XVolSize)*XVolSize
C But this is usually less efficient. On pipeline architectures it might however be better
         IF (Qx(I) .LT. 0.D0) THEN
            Qx(I) = Qx(I) + XVolSize
         ELSE
            IF (Qx(I) .GT. XVolSize) THEN
               Qx(I) = Qx(I) - XVolSize
            ENDIF
         ENDIF
         IF (Qy(I) .LT. 0.D0) THEN
            Qy(I) = Qy(I) + VolSize
         ELSE
            IF (Qy(I) .GT. VolSize) THEN
               Qy(I) = Qy(I) - VolSize
            ENDIF
         ENDIF
         IF (Qz(I) .LT. 0.D0) THEN
            Qz(I) = Qz(I) + VolSize
         ELSE
            IF (Qz(I) .GT. VolSize) THEN
               Qz(I) = Qz(I) - VolSize
            ENDIF
         ENDIF
      ENDDO

      CALL CalcPairList(Scale)
      CALL CalcForce (ForceX, ForceY, ForceZ, Virial)
      DO I=1,PartNum   
         Px(I) = Px(I) + 0.5D0*TimeStep*ForceX(I)
         Py(I) = Py(I) + 0.5D0*TimeStep*ForceY(I)
         Pz(I) = Pz(I) + 0.5D0*TimeStep*ForceZ(I)
      ENDDO
      CALL CalcTemp(TotKin)
      END





      SUBROUTINE CalcForce(ForceX, ForceY, ForceZ, Virial)

C *** Calculation of forces. We assume that the forces are a super- **
C *** position of central-symmetric forces between two particles.   ** 

C**********************Global Variables***************************
      include   "md.glb"
C*****************************************************************

      DOUBLE PRECISION ForceX(PartNum), ForceY(PartNum), 
     .       ForceZ(PartNum),
     .       R2, R4, R8, R14, ForceConst, Sx, Sy, Sz, 
     .       RMin2, Virial

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

           ForceConst = 24.d0*(2*R14-R8)
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






      SUBROUTINE CalcPotent(Potential)

C *** Calculation of total potential energy. We assume that the*******
C *** total potential energy can be written as a a superposition of **
C *** central-symmetric forces between two particles.  *************** 

C**********************Global Variables***************************
      include   "md.glb" 
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
        
          Potential = Potential + 4*(R12-R6)
        ENDIF
      ENDDO
      END

      SUBROUTINE CalcPairList(Scale)

C *** Calculation of total potential energy. We assume that the*******
C *** total potential energy can be written as a a superposition of **
C *** central-symmetric forces between two particles.  *************** 

C**********************Global Variables***************************
      include   "md.glb" 
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
        IF (.NOT.Scale) THEN
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
      include   "md.glb" 
C*****************************************************************

      DOUBLE PRECISION R2, Dx, Dy, Dz

      INTEGER I, J

      PairNum = 0
      DO I = 1,PartNum
        DO J=I+1,PartNum
C The following is not the most compact way of accounting for periodic 
C boundary conditions, but it is probably the most efficient.
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
           ENDIF
        ENDDO
      ENDDO
      IF (PairNum.GT.MaxLen) THEN
        print *, 'not enough memory allocated for pairlist'
      ENDIF
      END


      SUBROUTINE CalcTemp(KinEner)

C *** Calculation of total kinetic energy, i.e. the temperature ***
C *** Result is used for rescaling of velocities. *****************

C**********************Global Variables***************************
      include   "md.glb"
C*****************************************************************

      DOUBLE PRECISION KinEner

      INTEGER I

      KinEner = 0.D0
      DO I=1, PartNum
         KinEner = KinEner + Px(I)*Px(I) + Py(I)*Py(I)+
     .                       Pz(I)*Pz(I)
      ENDDO
      KinEner = KinEner*0.5D0
      END




      SUBROUTINE Rescale (KinEner)

C *** Rescaling velocities to adjust temperature. ****************

C**********************Global Variables***************************
C*                                                               *
      include  "md.glb"
C*                                                               *
C*                 Stored in file: "md.glb"                      *
C*  Px, Py, Pz, Qx, Qy, Qz : DOUBLE PRECISION arrays of size partnum       *
C*  PartNum, Temperature, Volume, InitTime, ScaleTime, TimeStep, *
C*  SimTime                                                      *
C*****************************************************************

      DOUBLE PRECISION KinEner, ScaleParam
      INTEGER I

C      return   

      KinEner = KinEner/(PartNum-1)
      ScaleParam = SQRT(Temperature*1.5D0/KinEner)

      DO I=1,PartNum
         Px(I) = Px(I)*ScaleParam
         Py(I) = Py(I)*ScaleParam
         Pz(I) = Pz(I)*ScaleParam
      ENDDO
      END






