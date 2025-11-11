      PROGRAM HarmOscDMC
C*****************************************************************
C* This program estimates the ground state energy of the         *
C* hydrogen atom using a diffusion quantum Monte Carlo           *
C* procedure with a guide function.                              *
C*                                                               *
C* Program written by Jos Thijssen	                         *
C* Summer 1999                                                   *
C* This program is described in section 12.2.2 of the book       *
C* "Computational Physics" by Jos Thijssen,                      *
C* Cambridge University Press 1999                               *
C*****************************************************************

      DOUBLE PRECISION Energy, DeltaT, Alpha

      INTEGER InitStep, TargetNum, StepNum

      CALL InitWalkers(InitStep, StepNum, Alpha, TargetNum, DeltaT)

      CALL Metropolis(Energy, InitStep, StepNum, Alpha, TargetNum,
     .                DeltaT)
      print *, Energy
      CLOSE(8)
      END



      SUBROUTINE InitWalkers(InitStep, StepNum, Alpha, TargetNum,
     .                       DeltaT)
C Initialise the parameters of the procedure and initialse the 
C random walkers to random initial positions. 

      include "globHydrDMC"

      INTEGER K, InitStep, StepNum, TargetNum, L

      DOUBLE PRECISION RealRand, DeltaT, Alpha

      print *, 'Give nr of steps and nr of equilibration steps'
      read *, StepNum, InitStep
      print *, 'Give target nr of Walkers'
      read *, TargetNum
      IF (TargetNum.GT.MaxWalkNum) THEN
        print *, 'This number should be smaller than', MaxWalkNum
        STOP
      END IF 
      print *, 'Give time step'
      read *, DeltaT
      print *, 'Give guide function parameter Alpha'
      read *, Alpha
      CALL InitRand(4537)
      DO K=1, TargetNum
        DO L=1, 3
          Walker(K, L) = RealRand()-0.5D0
        END DO
      END DO
      OPEN (8, File='ground_state')
      END



      DOUBLE PRECISION FUNCTION ELocal(R, Alpha)

      include "globHydrDMC"

      DOUBLE PRECISION R, Alpha

      ELocal = -1/R-0.5D0*Alpha*(Alpha-2.D0/R)
      END
      


      SUBROUTINE Metropolis (Energy, InitStep, StepNum, Alpha, 
     .                       TargetNum, DeltaT)
C The metropolis procedure. StepNum steps are carried out. In one step,
C all WalkNum walkers are moved to new positions. The first InitStep
C steps are used for equilibration. The energy is evaluated for the 
C variational parameter Alpha.  

      include "globHydrDMC"

      DOUBLE PRECISION Energy, Pos(3), NewPos(3),
     .                 GaussWidth, Rand1, Rand2, RandArr(6*MaxWalkNum), 
     .                 RealRand, DeltaT, Acc, 
     .                 R, NewR, EZero, Alpha, ELocal

      INTEGER K, Step, NAcc, I, Hist(0:100), NetStep, KK, L, 
     .        InitStep, StepNum, TargetNum, NewWN, NCopies, 
     .        WalkNum
      LOGICAL Alive(MaxWalkNum)
      
      Energy = 0.D0
      EZero = 0.D0
      NAcc = 0
      DeltaT = 0.05D0
      GaussWidth = SQRT(DeltaT)
      DO I=0,100
        Hist(I) = 0
      END DO
      DO I=1, MaxWalkNum
        Alive(I) = .TRUE.
      END DO
      WalkNum = TargetNum
      DO Step = 1, StepNum
        DO K=1, 3*WalkNum/2
C Store all random numbers needed for one time step in array RandArr...
C WalkNum should be EVEN!!!!
          CALL ExpRand(Rand1, Rand2)
          RandArr(2*K-1) = Rand1*GaussWidth
          RandArr(2*K)   = Rand2*GaussWidth
        END DO
        DO K=1, WalkNum
C Move walker nr. K and calculate the value R**2 for the 
C old and new value. Note that efficiency could be improved by 
C storing R**2 alongside the walker's position: in the present
C program, NewRSq or RSq is recalculated as RSq in the next step. 
          DO L=1, 3
            Pos(L)  = Walker(K, L)
          END DO
          R = SQRT(Pos(1)**2+Pos(2)**2+Pos(3)**2)
          DO L=1, 3
C F(L) = -2*Pos(L)*Alpha/R for trial function exp(-alpha r)....
            NewPos(L) = Pos(L) + RandArr(3*(K-1)+L) + 
     .                  DeltaT*(-Pos(L)*Alpha/R)
            Walker(K, L) = NewPos(L)
          END DO
          NewR = SQRT(NewPos(1)**2+NewPos(2)**2+NewPos(3)**2)
C There is no Metropolis correction built in. See the Helium atom prog for this
C Calculation of the acceptance rate
          Acc = EXP(DeltaT*(-0.5D0*(ELocal(NewR, Alpha)+
     .                              ELocal(R, Alpha))
     .                      +EZero))
C The acceptance rate is used to calculate the number
C of copies NCopies
          NCopies = INT(Acc+RealRand())
C Append NCopies-1 copies of Walker(K,:) to the end of the array Walker
          DO KK=1, NCopies-1
            WalkNum = WalkNum+1
            DO L=1, 3
              Walker(WalkNum, L) = NewPos(L)
            END DO
            Alive(WalkNum) = .TRUE.
          END DO
C Kill the walker if NCopies equals 0
          IF (NCopies.EQ.0) THEN
            Alive(K) = .FALSE.
          ELSE
            R = NewR
          END IF
c Update histogram...
          IF (Step.GT.InitStep) THEN
            IF (ABS(R).LT.4.D0) THEN
              Hist(NINT(R*10)) = Hist(NINT(R*10))+1
            END IF
          END IF
        ENDDO
C Compactify the array Walker, by removing all dead walkers
        NewWN = 0
        DO K=1, WalkNum
          IF (Alive(K)) THEN
            NewWN = NewWN+1
            DO L=1, 3
              Walker(NewWN, L) = Walker(K,L)
            END DO
            Alive(NewWN)=.TRUE.
          ENDIF
        ENDDO
        WalkNum = NewWN
c Estimate new value of EZero
        Ezero= EZero+0.1D0*LOG(DBLE(TargetNum)/WalkNum)
        IF (Step.GT.InitStep) THEN
          Energy = Energy+EZero
          IF (ABS(R).LT.2.D0) THEN
            Hist(NINT(R*10)) = Hist(NINT(R*10))+1
          END IF
        END IF        
      ENDDO
      NetStep = StepNum-InitStep
      Energy = Energy/(NetStep)
      DO I=0,40
        write (8, '(F10.5 F10.5)') I*0.1D0, 
     .                             Hist(I)/DBLE(WalkNum*NetStep)
      ENDDO
      END
      

