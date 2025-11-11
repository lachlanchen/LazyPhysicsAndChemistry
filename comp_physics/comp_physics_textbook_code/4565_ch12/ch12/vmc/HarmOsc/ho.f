      PROGRAM HarmOscVMC
C*****************************************************************
C* This program estimates the ground state energy of the         *
C* harmonic oscillator using a variational quantum Monta Carlo   *
C* procedure.                                                    *
C*                                                               *
C* Program written by Jos Thijssen	                         *
C* Spring 1999  february 2006                                    *
C* This program is described in section 12.2.2 of the book       *
C* "Computational Physics" by Jos Thijssen,                      *
C* Cambridge University Press 2007                               *
C*****************************************************************

      DOUBLE PRECISION Alpha, Energy

      INTEGER InitStep, WalkNum, StepNum

      CALL InitWalkers(Alpha, InitStep, StepNum, WalkNum)

      CALL Metropolis(Energy, Alpha, InitStep, StepNum, WalkNum)
      print *, Energy
      print *, 'analytical energy, var', 0.5*alpha+1/(8*alpha), 
     .                             (1-4*alpha**2)**2/(32*alpha**2)
      CLOSE(11)
      CLOSE(12)
      CLOSE(8)
      END



      SUBROUTINE InitWalkers(Alpha, InitStep, StepNum, WalkNum)
C Initialise the parameters of the procedure and initialse the 
C random walkers to random initial positions. 

      include "globHO"

      INTEGER K, InitStep, StepNum, WalkNum

      DOUBLE PRECISION RealRand, Alpha

      print *, 'Give nr of steps and nr of equilibration steps'
      read *, StepNum, InitStep
      print *, 'Give nr of Walkers'
      read *, WalkNum
      IF (WalkNum.GT.MaxWalkNum) THEN
        print *, 'This number should be smaller than', MaxWalkNum
        STOP
      END IF 
      print *, 'Give exponential parameter alpha'
      read *, Alpha
      CALL InitRand(4537)
      DO K=1, WalkNum
        Walker(K) = RealRand()-0.5D0
      ENDDO
      OPEN (11, File='energy.dat')
      OPEN (12, File='variance.dat')
      OPEN (8, File='ground_state.dat')
      END


      DOUBLE PRECISION FUNCTION TransProb(X, NewX, Alpha)
C Calculate the transition probability for a walker with coordinate 
C X stepping to NewX, with variational parameter Alpha

      include "globHO"

      DOUBLE PRECISION Alpha, X, NewX

      TransProb = EXP(-2*Alpha*(NewX**2-X**2))

      END




      SUBROUTINE CalcLocal (ELocal, X, Alpha)
C Calculate the local energy for a walker at position X

      include "globHO"

      DOUBLE PRECISION ELocal, Alpha, X

      ELocal = Alpha + X**2*(0.5D0-2.D0*Alpha**2)
      END
      


      SUBROUTINE Metropolis (Energy, Alpha, InitStep, StepNum, WalkNum)
C The metropolis procedure. StepNum steps are carried out. In one step,
C all WalkNum walkers are moved to new positions. The first InitStep
C steps are used for equilibration. The energy is evaluated for the 
C variational parameter Alpha.  

      include "globHO"

      DOUBLE PRECISION Energy, Alpha, X, NewX, TransProb, ELocal,
     .                 GaussWidth, Rand1, Rand2, RandArr(MaxWalkNum), 
     .                 RealRand, EnerSq

      INTEGER K, Step, NAcc, I, Hist(-100:100), NetStep,
     .        InitStep, StepNum, WalkNum
      
      Energy = 0.D0
      EnerSq = 0.D0
      NAcc = 0
      GaussWidth = 1.D0
      DO I=-100,100
        Hist(I) = 0
      END DO
      DO Step = 1, StepNum
        DO K=1, WalkNum/2
          CALL ExpRand(Rand1, Rand2)
          RandArr(2*K-1) = Rand1
          RandArr(2*K)   = Rand2
        END DO
        Energy = 0.D0
        EnerSq = 0.D0
        DO K=1, WalkNum
          X    = Walker(K)
          NewX = Walker(K) + RandArr(K)*GaussWidth
C Again, the expression for the transition probability is simple as a 
C result of the Gaussian form of the trial function
          IF (RealRand().LT.TransProb(X, NewX, Alpha)) THEN
            NAcc = NAcc+1
            Walker(K) = NewX
            X = NewX
          END IF
          IF (Step.GT.InitStep) THEN
            CALL CalcLocal(ELocal, X, Alpha)
            Energy = Energy + ELocal
            EnerSq = EnerSq + ELocal*ELocal
            IF (ABS(X).LT.2.D0) THEN
              Hist(NINT(X*5)) = Hist(NINT(X*5))+1
            END IF
          END IF
        ENDDO
C Dynamic adaptation of the width of the Gaussian distribution
        IF (MOD(Step,100) .EQ. 0) THEN
          GaussWidth = GaussWidth*NAcc/(50.D0*WalkNum)
          NAcc = 0
        END IF
        IF (Step.GT.InitStep) THEN
          Energy = Energy/WalkNum
          WRITE (11, *) Energy
          WRITE (12, *) EnerSq/WalkNum-Energy**2
          Energy = 0.d0
          EnerSq = 0.d0
        END IF
      ENDDO
      NetStep = StepNum-InitStep
      Energy = Energy/(NetStep*WalkNum)
      DO I=-10,10
        write (8, '(F10.5 F10.5)') I*0.2D0, 
     .                             Hist(I)/DBLE(WalkNum*NetStep)
      ENDDO
      END
      

