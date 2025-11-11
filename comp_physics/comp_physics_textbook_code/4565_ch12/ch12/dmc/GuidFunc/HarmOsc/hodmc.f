      PROGRAM HarmOscDMC
C*****************************************************************
C* This program estimates the ground state energy of the         *
C* harmonic oscillator using a diffusion quantum Monte Carlo     *
C* procedure with a guide function.                              *                                                  C*                                                               *
C* Program written by Jos Thijssen	                         *
C* Summer 1999                                                   *
C* This program is described in section 12.3.3 of the book       *
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

      include "globHODMC"

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



      DOUBLE PRECISION FUNCTION ELocal(RSq, Alpha)

      include "globHODMC"

      DOUBLE PRECISION RSq, Alpha

      ELocal = 3*Alpha+RSq*(0.5D0-2*Alpha**2)
      END
      


      SUBROUTINE Metropolis (Energy, InitStep, StepNum, Alpha, 
     .                       TargetNum, DeltaT)
C The metropolis procedure. StepNum steps are carried out. In one step,
C all WalkNum walkers are moved to new positions. The first InitStep
C steps are used for equilibration. The energy is evaluated for the 
C variational parameter Alpha.  
C The efficiency of this program can be improved considerably by
C storing extra quantities for each walker, such as the local energy
C and the drift term F. 

      include "globHODMC"

      DOUBLE PRECISION Energy, Pos(3), NewPos(3),
     .                 GaussWidth, Rand1, Rand2, RandArr(6*MaxWalkNum), 
     .                 RealRand, DeltaT, Acc, Weight, 
     .                 RSq, NewRSq, R, EZero, Alpha, ELocal,
     .                 OldF(3), NewF(3), Q

      INTEGER K, Step, NAcc, I, Hist(0:100), NetStep, KK, L, 
     .        InitStep, StepNum, TargetNum, NewWN, NCopies, 
     .        WalkNum
      LOGICAL Alive(MaxWalkNum)
      
      Energy = 0.D0
C The starting value should be reasonable, otherwise strong 
C oscillations arise in the population size....
      EZero = 1.5D0
      NAcc = 0
      GaussWidth = SQRT(DeltaT)
      DO I=0,100
        Hist(I) = 0
      END DO
      DO I=1, TargetNum
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
          RSq = 0.D0
          NewRSq = 0.D0
          DO L=1, 3
            Pos(L)  = Walker(K, L)
            NewPos(L) = Pos(L) + RandArr(3*(K-1)+L) + 
     .                  DeltaT*(-2*Pos(L)*Alpha)
            RSq = RSq + Pos(L)**2
            NewRSq = NewRSq + NewPos(L)**2
c             Walker(K, L) = NewPos(L)
            OldF(L) = -4*Alpha*Pos(L)
            NewF(L) = -4*Alpha*NewPos(L)
          END DO
C Metropolis step.....
          Q = 0.D0
          DO L=1, 3
            Q = Q + 0.5D0*(OldF(L)+NewF(L))*(0.25D0*DeltaT*
     .                (OldF(L)-NewF(L))-NewPos(L)+Pos(L))
          ENDDO
          Q = Q-2*Alpha*(NewRSq-RSq)
c          Q = 0
          Acc = min(1.d0, EXP(Q))
          IF (RealRand().LT.Acc) THEN
            DO L=1, 3
              Walker(K, L) = NewPos(L) 
            ENDDO
            Weight = EXP(DeltaT*(-0.5D0*(ELocal(NewRSq, Alpha)+
     .                              ELocal(RSq, Alpha))
     .                      +EZero))
          ELSE
            Weight = 1.d0
C Weight=1.d0 ensures that nothing changes....           
          ENDIF
C The acceptance rate is used to calculate the number
C of copies NCopies
          NCopies = INT(Weight+RealRand())
C Append NCopies-1 copies of Walker(K,:) to the end of the array Walker
          DO KK=1, NCopies-1
            IF (WalkNum.LT.MaxWalkNum) THEN 
              WalkNum = WalkNum+1
              DO L=1, 3
                Walker(WalkNum, L) = NewPos(L)
              END DO
              Alive(WalkNum) = .TRUE.
            ELSE
              IF (Step.GT.InitStep) THEN
                print *, 'too many walkers..'
                stop
              END IF
            END IF
          END DO
C Kill the walker if NCopies equals 0
            IF (NCopies.EQ.0) THEN
              Alive(K) = .FALSE.
            ELSE
              RSq = NewRSq
            END IF
c Update histogram...
          IF (Step.GT.InitStep) THEN
            R = SQRT(RSq)
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
        IF (MOD(Step,100).EQ.0) print *, NewWN
c Estimate new value of EZero
        Ezero = EZero+0.05D0*LOG(DBLE(TargetNum)/WalkNum)
        IF (Step.GT.InitStep) THEN
          Energy = Energy+EZero
          R = SQRT(RSq)
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
      

